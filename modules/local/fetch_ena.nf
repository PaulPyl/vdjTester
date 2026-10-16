// Download one 10x library from ENA (one or more runs = lanes) and rename the
// files to Cell Ranger's convention: <fastq_id>_S1_L00<lane>_R<read>_001.fastq.gz.
//
// ENA serves SRA's regenerated _1/_2 files without index reads, and the order is
// not guaranteed (in GSE266035 `_1` holds the cDNA read). Read length does not
// tell them apart either: VDJ and GEX libraries are often sequenced 150 x 150.
// In 10x 5' chemistry read 1 is barcode (16) + UMI (10) + template-switch oligo
// (TTTCTTATATGGG from position 27). The exact oligo is often degraded (cll_5
// timepoint 3 BCR: 7 % exact matches), so the test looks for its core TCTTATAT in
// bases 23-42: 43-70 % of read 1 vs <= 0.1 % of read 2 across GSE311112 runs.
//
// File URLs and MD5s come from the ENA portal API. When the API is down (it
// returned HTTP 500 and timeouts on 2026-10-16), the URL is built from ENA's
// fixed layout and the MD5 check is skipped, which is logged in fetch.log.
//
// Raw files are downloaded into <store_dir>/fastq/.downloads rather than the work
// dir, so a failed or killed attempt resumes where it stopped (a GEX library is
// 6-13 GB, and ENA served ~1 MB/s per file on 2026-10-16). A cached file whose
// MD5 already matches is not fetched again; one that mismatches is deleted.
process FETCH_ENA {
    tag "${fastq_id}"
    label 'process_single'
    storeDir "${params.store_dir}/fastq"

    input:
    tuple val(sample_id), val(fastq_id), val(runs)

    output:
    tuple val(sample_id), val(fastq_id), path(fastq_id)

    script:
    """
    fetch() {  # url, output file; resumes after dropped connections
        # No curl --retry here: after an HTTP 500 mid-download, curl's own retry
        # resumed from a stale offset and the file failed its MD5 (SRR36170693,
        # 2026-10-16). A fresh curl per attempt takes the offset from the file size.
        for attempt in \$(seq 1 30); do
            curl -fsSL -C - -o "\$2" "\$1" && return 0
            echo "download of \$2 interrupted (attempt \$attempt), resuming" >&2
            sleep 30
        done
        return 1
    }
    md5_of() { if command -v md5sum >/dev/null; then md5sum "\$1" | cut -d' ' -f1; else md5 -q "\$1"; fi; }
    tso_pct() {  # % of the first 2000 reads with the TSO core just after barcode + UMI
        set +o pipefail
        gzip -dc "\$1" | head -n 8000 | awk 'NR % 4 == 2 { n++; if (index(substr(\$0, 23, 20), "TCTTATAT")) t++ } END { print int(100 * t / n) }'
    }
    ena_dir() {  # ENA's layout: vol1/fastq/<first 6>/<00x|0xx|xxx>/<acc>
        acc=\$1; digits=\${acc#???}; n=\${#digits}
        case \$n in
            6) sub="" ;;
            7) sub="00\${acc: -1}/" ;;
            8) sub="0\${acc: -2}/" ;;
            *) sub="\${acc: -3}/" ;;
        esac
        echo "ftp.sra.ebi.ac.uk/vol1/fastq/\${acc:0:6}/\${sub}\${acc}"
    }

    cache=${params.store_dir}/fastq/.downloads
    out=${fastq_id}.partial
    mkdir -p "\$out" "\$cache"
    lane=0
    for run in ${runs.join(' ')}; do
        lane=\$((lane + 1))
        report=\$(curl -fsSL --max-time 60 --retry 3 --retry-all-errors \\
            "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=\$run&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" \\
            | tail -n 1 || true)
        urls=\$(echo "\$report" | cut -f2 | tr ';' ' ')
        md5s=\$(echo "\$report" | cut -f3 | tr ';' ' ')
        if [ -z "\$urls" ] || ! echo "\$urls" | grep -q "\$run"; then
            echo "\$run: ENA API unavailable, using the fixed path layout without MD5 check" >> "\$out/fetch.log"
            d=\$(ena_dir "\$run")
            urls="\$d/\${run}_1.fastq.gz \$d/\${run}_2.fastq.gz"
            md5s=""
        fi
        set -- \$md5s
        i=0
        for url in \$urls; do
            i=\$((i + 1))
            f="\$cache/\$run.\$i.fastq.gz"
            expected=""
            [ \$# -ge \$i ] && expected=\$(eval echo "\\\${\$i}")
            if [ -n "\$expected" ] && [ -s "\$f" ] && [ "\$(md5_of "\$f")" = "\$expected" ]; then
                echo "\$run file \$i: already in the download cache" >> "\$out/fetch.log"
                continue
            fi
            fetch "https://\$url" "\$f"
            if [ -n "\$expected" ] && [ "\$(md5_of "\$f")" != "\$expected" ]; then
                echo "\$run file \$i: md5 mismatch, downloading it again from the start" >> "\$out/fetch.log"
                rm -f "\$f"
                fetch "https://\$url" "\$f"
                actual=\$(md5_of "\$f")
                if [ "\$actual" != "\$expected" ]; then
                    echo "md5 mismatch for \$run file \$i: \$actual vs \$expected" >&2
                    rm -f "\$f"
                    exit 1
                fi
            fi
        done
        if [ "\$i" -ne 2 ]; then
            echo "\$run: expected 2 FASTQ files from ENA, got \$i" >&2
            exit 1
        fi

        tso1=\$(tso_pct "\$cache/\$run.1.fastq.gz"); tso2=\$(tso_pct "\$cache/\$run.2.fastq.gz")
        if [ "\$tso1" -ge 20 ] && [ "\$tso2" -lt 5 ]; then
            r1=\$cache/\$run.1.fastq.gz; r2=\$cache/\$run.2.fastq.gz
        elif [ "\$tso2" -ge 20 ] && [ "\$tso1" -lt 5 ]; then
            r1=\$cache/\$run.2.fastq.gz; r2=\$cache/\$run.1.fastq.gz
        else
            echo "\$run: cannot tell R1 from R2 (TSO core in bases 23-42: \$tso1 % / \$tso2 %)" >&2
            exit 1
        fi
        echo "\$run (lane \$lane): TSO core in bases 23-42 of \$tso1 % / \$tso2 % of reads, R1 = \$(basename "\$r1")" >> "\$out/fetch.log"
        # Moved out of the cache only once every lane is in, so a failure in a
        # later lane does not lose the earlier ones.
        echo "\$r1 \$out/${fastq_id}_S1_L00\${lane}_R1_001.fastq.gz" >> moves.txt
        echo "\$r2 \$out/${fastq_id}_S1_L00\${lane}_R2_001.fastq.gz" >> moves.txt
    done
    while read -r src dst; do mv "\$src" "\$dst"; done < moves.txt
    mv "\$out" ${fastq_id}
    """

    stub:
    """
    mkdir -p ${fastq_id}
    lane=0
    for run in ${runs.join(' ')}; do
        lane=\$((lane + 1))
        touch ${fastq_id}/${fastq_id}_S1_L00\${lane}_R1_001.fastq.gz ${fastq_id}/${fastq_id}_S1_L00\${lane}_R2_001.fastq.gz
    done
    """
}
