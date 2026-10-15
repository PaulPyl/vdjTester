// Download a tarball once and keep it unpacked in params.store_dir across runs.
// The top-level directory of the archive is stripped, so the output is always
// a directory called `name`.
process FETCH_TAR {
    tag "${name}"
    label 'process_single'
    storeDir "${params.store_dir}"

    input:
    tuple val(key), val(name), val(url)

    output:
    tuple val(key), path(name)

    script:
    def gz = url.endsWith('.gz') ? 'z' : ''
    """
    set -o pipefail
    mkdir -p ${name}.partial
    curl -fsSL --retry 3 "${url}" | tar -x${gz}f - -C ${name}.partial --strip-components=1
    mv ${name}.partial ${name}
    """

    stub:
    """
    mkdir -p ${name}
    """
}
