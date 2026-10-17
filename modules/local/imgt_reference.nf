// The IMGT reference option (--reference_set imgt): every method on IMGT germlines,
// each through its own published route, with nothing rebuilt or filtered here.
//   Cell Ranger   fetch-imgt (shipped with Cell Ranger) + cellranger mkvdjref, the
//                 route 10x documents for IMGT references; IMGT's current release
//   MiXCR         MiLaboratories' IMGT library (github.com/repseqio/library-imgt),
//                 params.mixcr_imgt_library; fetched with FETCH_FILES in main.nf
//   Immcantation  the IMGT germlines in its container, as in the default run
//   TRUST4        the IMGT file of its release, as in the default run
// The IMGT releases therefore differ between methods (see docs/PLAN.md); what this
// option removes is the non-IMGT reference of Cell Ranger (Ensembl 94) and MiXCR's
// own library.
//
// IMGT serves only its current release, so the fetch-imgt download is kept in
// <store_dir>/imgt and reused; delete that directory to refresh it.

// fetch-imgt needs Cell Ranger's Python; IMGT's query server answers 504 now and
// then, so it is retried.
process FETCH_IMGT {
    label 'process_single'
    storeDir "${params.store_dir}/imgt"

    output:
    path 'imgt_fetch'

    script:
    def cr = params.cellranger_dir
    """
    mkdir -p imgt_fetch && cd imgt_fetch
    for attempt in 1 2 3 4 5; do
        PYTHONPATH=${cr}/lib/python ${cr}/external/anaconda/bin/python ${cr}/lib/bin/fetch-imgt --genome imgt > fetch-imgt.log 2>&1 \\
            && [ -s imgt-mkvdjref-input.fasta ] && break
        echo "fetch-imgt attempt \$attempt failed: \$(tail -n 1 fetch-imgt.log)" >&2
        sleep 60
    done
    [ -s imgt-mkvdjref-input.fasta ] || { echo "fetch-imgt failed 5 times" >&2; exit 1; }
    rm -f *.html
    date -u +%Y-%m-%d > downloaded.txt
    """

    stub:
    """
    mkdir -p imgt_fetch && touch imgt_fetch/imgt-mkvdjref-input.fasta imgt_fetch/downloaded.txt
    """
}

process MKVDJREF_IMGT {
    label 'process_single'
    storeDir "${params.store_dir}/imgt"

    input:
    path fetched

    output:
    path 'vdj_imgt'

    script:
    def cr_path = params.cellranger_dir ? "export PATH=${params.cellranger_dir}:\$PATH" : ''
    """
    ${cr_path}
    cellranger mkvdjref --genome=vdj_imgt --seqs=${fetched}/imgt-mkvdjref-input.fasta
    cp ${fetched}/downloaded.txt vdj_imgt/imgt_downloaded.txt
    rm -rf mkvdjref_vdj_imgt
    """

    stub:
    """
    mkdir -p vdj_imgt/fasta && touch vdj_imgt/fasta/regions.fa vdj_imgt/reference.json vdj_imgt/imgt_downloaded.txt
    """
}
