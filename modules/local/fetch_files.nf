// Download a few small files once into params.store_dir/<name>/, e.g. the TRUST4
// germline references, which are not part of its container.
process FETCH_FILES {
    tag "${name}"
    label 'process_single'
    storeDir "${params.store_dir}"

    input:
    tuple val(name), val(files)   // [[file name, url], ...]

    output:
    path name

    script:
    // Explicit names: curl -O would keep URL escapes (human_IMGT%2BC.fa).
    def fetches = files.collect { f -> "curl -fsSL --retry 5 --retry-all-errors -o '${f[0]}' '${f[1]}'" }.join('\n    ')
    """
    mkdir -p ${name}.partial
    cd ${name}.partial
    ${fetches}
    cd ..
    mv ${name}.partial ${name}
    """

    stub:
    """
    mkdir -p ${name}
    """
}
