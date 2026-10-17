// The unified germline reference (--reference_set imgt): every tool's V, D and J
// germlines from one IMGT snapshot, instead of each tool's own default (10x's
// Ensembl-94 reference, Immcantation's IMGT 2022-12-16, MiXCR's built-in library,
// TRUST4's IMGT file). See bin/prep_imgt_refs.py for which alleles are kept.
//
// IMGT only serves its current release, so the snapshot is downloaded once into
// <store_dir>/imgt and reused; delete that directory to refresh it. snapshot.tsv
// records the download date and the bulk file's Last-Modified date.

// One IMGT snapshot: the bulk GENE-DB file, and Cell Ranger's fetch-imgt queries
// (leader + V, D, J, C), downloaded together. fetch-imgt needs Cell Ranger's Python;
// IMGT's query server answers 504 now and then, so it is retried.
process FETCH_IMGT {
    label 'process_single'
    storeDir "${params.store_dir}/imgt"

    output:
    path 'imgt_snapshot'

    script:
    // Cell Ranger's bundled Python does the downloads too: the container that runs
    // Cell Ranger on the Mac (plain Ubuntu) has no curl.
    def cr = params.cellranger_dir
    def py = "PYTHONPATH=${cr}/lib/python ${cr}/external/anaconda/bin/python"
    """
    mkdir -p imgt_snapshot && cd imgt_snapshot
    url='${params.imgt_bulk_url}'
    for attempt in 1 2 3 4 5; do
        modified=\$(${py} - "\$url" <<'PY'
import shutil, ssl, sys, urllib.request
import certifi
ctx = ssl.create_default_context(cafile=certifi.where())
with urllib.request.urlopen(sys.argv[1], context=ctx, timeout=300) as r, open("bulk.fa", "wb") as out:
    shutil.copyfileobj(r, out)
    print(r.headers.get("Last-Modified", ""))
PY
        ) && [ -s bulk.fa ] && break
        echo "IMGT bulk download attempt \$attempt failed" >&2
        sleep 60
    done
    [ -s bulk.fa ] || { echo "IMGT bulk download failed 5 times" >&2; exit 1; }

    for attempt in 1 2 3 4 5; do
        ${py} ${cr}/lib/bin/fetch-imgt --genome imgt > fetch-imgt.log 2>&1 \\
            && [ -s imgt-mkvdjref-input.fasta ] && break
        echo "fetch-imgt attempt \$attempt failed: \$(tail -n 1 fetch-imgt.log)" >&2
        sleep 60
    done
    [ -s imgt-mkvdjref-input.fasta ] || { echo "fetch-imgt failed 5 times" >&2; exit 1; }
    mv imgt-mkvdjref-input.fasta mkvdjref-input.fasta
    mv imgt-imgt-raw.fasta fetch-imgt-raw.fasta
    rm -f *.html

    printf 'key\\tvalue\\n' > snapshot.tsv
    printf 'downloaded\\t%s\\n' "\$(date -u +%Y-%m-%d)" >> snapshot.tsv
    printf 'imgt_bulk_last_modified\\t%s\\n' "\$modified" >> snapshot.tsv
    printf 'imgt_bulk_url\\t%s\\n' "\$url" >> snapshot.tsv
    """

    stub:
    """
    mkdir -p imgt_snapshot
    touch imgt_snapshot/bulk.fa imgt_snapshot/mkvdjref-input.fasta imgt_snapshot/snapshot.tsv
    """
}

// Per-tool inputs from the snapshot, restricted to the alleles both parts have.
process PREP_IMGT_REFS {
    label 'process_single'
    storeDir "${params.store_dir}/imgt"

    input:
    path snapshot
    path trust4_default

    output:
    path 'imgt_refs'

    script:
    """
    mkdir -p imgt_refs
    prep_imgt_refs.py \\
        --bulk ${snapshot}/bulk.fa \\
        --cellranger-input ${snapshot}/mkvdjref-input.fasta \\
        --trust4-default ${trust4_default}/human_IMGT+C.fa \\
        --outdir imgt_refs > imgt_refs/summary.txt
    cp ${snapshot}/snapshot.tsv imgt_refs/
    # TRUST4 also needs its genome-coordinate candidate FASTA (read recruitment only).
    cp ${trust4_default}/hg38_bcrtcr.fa imgt_refs/trust4/
    """

    stub:
    """
    mkdir -p imgt_refs/cellranger imgt_refs/germlines/human/vdj imgt_refs/mixcr imgt_refs/trust4
    touch imgt_refs/cellranger/seqs.fasta imgt_refs/trust4/human_IMGT+C.fa imgt_refs/trust4/hg38_bcrtcr.fa imgt_refs/snapshot.tsv imgt_refs/summary.txt
    """
}

process MKVDJREF_IMGT {
    label 'process_single'
    storeDir "${params.store_dir}/imgt"

    input:
    path refs

    output:
    path 'vdj_imgt'

    script:
    def cr_path = params.cellranger_dir ? "export PATH=${params.cellranger_dir}:\$PATH" : ''
    """
    ${cr_path}
    cellranger mkvdjref --genome=vdj_imgt --seqs=${refs}/cellranger/seqs.fasta
    rm -rf mkvdjref_vdj_imgt
    """

    stub:
    """
    mkdir -p vdj_imgt/fasta && touch vdj_imgt/fasta/regions.fa vdj_imgt/reference.json
    """
}

// IgBLAST database for changeo-10x (-b), next to the gapped germlines (-r).
// - internal_data (FWR/CDR positions) comes from the container: IgBLAST takes domain
//   annotation from its own internal database, whatever the germline set.
// - C-region databases come from the container: constant genes are not compared.
// - The J aux file (frame, CDR3 end) is IgBLAST's own line for every J allele it
//   knows, else the line prep_imgt_refs.py derived from the sequence. Without a line
//   IgBLAST reports no junction: on cll_1 every lambda contig (IGLJ2A/2B, newer IMGT
//   names) was dropped. The derived lines match IgBLAST's for 114 of 126 alleles; the
//   rest are non-canonical J motifs IgBLAST curated by hand.
process IGBLAST_IMGT {
    label 'process_single'
    storeDir "${params.store_dir}/imgt"

    input:
    path refs

    output:
    path 'igblast_imgt'

    script:
    """
    mkdir -p igblast_imgt
    imgt2igblast.sh -i ${refs}/germlines -o igblast_imgt > igblast_imgt/imgt2igblast.log 2>&1
    cp /usr/local/share/igblast/database/imgt_human_ig_c.* /usr/local/share/igblast/database/imgt_human_tr_c.* igblast_imgt/database/
    cp -r /usr/local/share/igblast/internal_data igblast_imgt/
    mkdir -p igblast_imgt/optional_file
    awk 'NR == FNR { if (NF && \$1 !~ /^#/) known[\$1] = \$0; next }
         !NF || /^#/ { print; next }
         { print (\$1 in known) ? known[\$1] : \$0 }' \\
        /usr/local/share/igblast/optional_file/human_gl.aux ${refs}/igblast/human_gl.aux \\
        > igblast_imgt/optional_file/human_gl.aux
    cp -r ${refs}/germlines igblast_imgt/
    """

    stub:
    """
    mkdir -p igblast_imgt/germlines/human/vdj igblast_imgt/database
    """
}

// MiXCR library: one buildLibrary per chain (anchor points inferred by alignment to
// the built-in library), merged. TRAV/DV genes rearrange in both TRA and TRD, as in
// MiXCR's own library; C genes are MiXCR's built-in ones.
process MIXCR_LIBRARY_IMGT {
    label 'process_single'
    storeDir "${params.store_dir}/imgt"
    secret 'MI_LICENSE'

    input:
    path refs

    output:
    path 'mixcr_imgt.json.gz'

    script:
    """
    export HOME=\$PWD TMPDIR=\$PWD XDG_CACHE_HOME=\$PWD/.cache
    cp ${refs}/mixcr/*.fasta .
    libs=""
    for c in IGH IGK IGL TRA TRB TRD TRG; do
        lc=\$(echo \$c | tr 'A-Z' 'a-z')
        d=""; [ -s \${c}D.fasta ] && d="--d-genes-from-fasta \${c}D.fasta"
        mixcr buildLibrary --force-overwrite --species hsa --chain \$c --taxon-id 9606 \\
            --v-genes-from-fasta \${c}V.fasta --v-gene-feature VRegion \\
            --j-genes-from-fasta \${c}J.fasta \$d \\
            --c-genes-from-species hsa \$lc.json.gz > build_\$c.log 2>&1
        libs="\$libs \$lc.json.gz"
    done
    if command -v python3 > /dev/null; then python3 - <<'PY'
import gzip, json
lib = json.load(gzip.open("tra.json.gz"))
for entry in lib:
    for gene in entry["genes"]:
        if "/DV" in gene["name"]:
            gene["chains"] = ["TRA", "TRD"]
json.dump(lib, gzip.open("tra.json.gz", "wt"))
PY
    else
        echo "WARNING: no python3, TRAV/DV genes stay TRA-only in the MiXCR library" >&2
    fi
    mixcr mergeLibrary --force-overwrite \$libs mixcr_imgt.json.gz
    """

    stub:
    """
    touch mixcr_imgt.json.gz
    """
}
