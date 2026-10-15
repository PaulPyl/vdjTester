#!/bin/bash
# Launch vdjTester on Dardel. Run from a login node inside tmux/screen:
#   ./run_dardel.sh                 # full run
#   ./run_dardel.sh -stub           # wiring check; extra args go to nextflow
#
# Prerequisites (once per host):
#   * Cell Ranger unpacked at params.cellranger_dir (conf/dardel.config)
#   * nextflow secrets set MI_LICENSE <key>   (stored in ~/.nextflow, so NXF_HOME
#     is deliberately left at its default here)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

module load PDC/26.03 apptainer/1.5.1-cpeGNU-26.03
module load nextflow/26.04.3

# Keep Apptainer's own caches in the project, never in $HOME. nextflow.config can
# only set the image cache; the blob cache and conversion tmpdir belong to the
# apptainer process. The tmpdir needs room for a whole uncompressed image.
export APPTAINER_CACHEDIR="${PWD}/container_cache/apptainer"
export APPTAINER_TMPDIR="${PWD}/container_cache/tmp"
export NXF_APPTAINER_CACHEDIR="${PWD}/container_cache/nextflow"
export NXF_TEMP="${PWD}/container_cache/tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$NXF_APPTAINER_CACHEDIR"

if ! nextflow secrets list 2>/dev/null | grep -qx MI_LICENSE; then
    echo "ERROR: Nextflow secret MI_LICENSE is not set (nextflow secrets set MI_LICENSE <key>)" >&2
    exit 1
fi

# Stub outputs are empty placeholders. Keep them out of the real store_dir, or a
# later real run would take them for finished downloads.
extra=()
for arg in "$@"; do
    if [ "$arg" = "-stub" ] || [ "$arg" = "-stub-run" ]; then
        extra=(--store_dir "${PWD}/results_stub/store" --outdir "${PWD}/results_stub")
    fi
done

# -ansi-log false: plain lines, since the output is piped into a log file.
nextflow run . -profile dardel -resume -ansi-log false "${extra[@]}" "$@" 2>&1 | tee "run_$(date +%Y%m%d_%H%M%S).log"
