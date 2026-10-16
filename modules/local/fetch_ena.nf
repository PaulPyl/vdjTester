// Download one 10x library from ENA (one or more runs = lanes) and rename the
// files to Cell Ranger's convention: <fastq_id>_S1_L00<lane>_R<read>_001.fastq.gz.
//
// ENA serves SRA's regenerated _1/_2 files without index reads, and the order is
// not guaranteed (in GSE266035 `_1` holds the cDNA read). Read length does not
// tell them apart either: VDJ and GEX libraries are often sequenced 150 x 150.
// In 10x 5' chemistry read 1 is barcode (16) + UMI (10) + template-switch oligo,
// so the read with TTTCTTATATGGG at position 27 is R1 (~60 % of reads in
// GSE311112 vs 0 % in the other file).
//
// File URLs and MD5s come from the ENA portal API. When the API is down (it
// returned HTTP 500 and timeouts on 2026-10-16), the URL is built from ENA's
// fixed layout and the MD5 check is skipped, which is logged in fetch.log.
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
        for attempt in 1 2 3 4 5 6 7 8 9 10; do
            curl -fSL --retry 5 --retry-all-errors -C - -o "\$2" "\$1" && return 0
            echo "download of \$2 interrupted (attempt \$attempt), resuming" >&2
            sleep 15
        done
        return 1
    }
    md5_of() { if command -v md5sum >/dev/null; then md5sum "\$1" | cut -d' ' -f1; else md5 -q "\$1"; fi; }
    tso_pct() {  # % of the first 2000 reads with the 5' TSO right after barcode + UMI
        set +o pipefail
        gzip -dc "\$1" | head -n 8000 | awk 'NR % 4 == 2 { n++; if (substr(\$0, 27, 13) == "TTTCTTATATGGG") t++ } END { print int(100 * t / n) }'
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

    out=${fastq_id}.partial
    mkdir -p "\$out"
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
            fetch "https://\$url" "\$run.\$i.fastq.gz"
            if [ \$# -ge \$i ]; then
                expected=\$(eval echo "\\\${\$i}")
                actual=\$(md5_of "\$run.\$i.fastq.gz")
                if [ "\$actual" != "\$expected" ]; then
                    echo "md5 mismatch for \$run file \$i: \$actual vs \$expected" >&2
                    exit 1
                fi
            fi
        done
        if [ "\$i" -ne 2 ]; then
            echo "\$run: expected 2 FASTQ files from ENA, got \$i" >&2
            exit 1
        fi

        tso1=\$(tso_pct "\$run.1.fastq.gz"); tso2=\$(tso_pct "\$run.2.fastq.gz")
        if [ "\$tso1" -ge 20 ] && [ "\$tso2" -lt 5 ]; then
            r1=\$run.1.fastq.gz; r2=\$run.2.fastq.gz
        elif [ "\$tso2" -ge 20 ] && [ "\$tso1" -lt 5 ]; then
            r1=\$run.2.fastq.gz; r2=\$run.1.fastq.gz
        else
            echo "\$run: cannot tell R1 from R2 (TSO at position 27: \$tso1 % / \$tso2 %)" >&2
            exit 1
        fi
        echo "\$run (lane \$lane): TSO at position 27 in \$tso1 % / \$tso2 % of reads, R1 = \$r1" >> "\$out/fetch.log"
        mv "\$r1" "\$out/${fastq_id}_S1_L00\${lane}_R1_001.fastq.gz"
        mv "\$r2" "\$out/${fastq_id}_S1_L00\${lane}_R2_001.fastq.gz"
    done
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
