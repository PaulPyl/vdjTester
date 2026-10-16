// Download a tarball once and keep it unpacked in params.store_dir across runs.
// The top-level directory of the archive is stripped, so the output is always
// a directory called `name`.
//
// The archive goes to a file first, not straight into tar: 10x's CDN reset the
// connection partway through the 11.5 GB GEX reference (2026-10-16), and a
// stream cannot resume. `curl -C -` continues from where the file stopped.
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
    archive=${name}.download
    for attempt in 1 2 3 4 5 6 7 8 9 10; do
        if curl -fSL --retry 5 --retry-all-errors -C - -o "\$archive" "${url}"; then
            break
        fi
        if [ "\$attempt" -eq 10 ]; then
            echo "download failed after \$attempt attempts" >&2
            exit 1
        fi
        echo "download interrupted (attempt \$attempt), resuming" >&2
        sleep 15
    done

    mkdir -p ${name}.partial
    tar -x${gz}f "\$archive" -C ${name}.partial --strip-components=1
    rm -f "\$archive"
    mv ${name}.partial ${name}
    """

    stub:
    """
    mkdir -p ${name}
    """
}
