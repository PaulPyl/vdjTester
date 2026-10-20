#!/bin/bash
#SBATCH -A naiss2026-4-384-cpu
#SBATCH -p cpu
#SBATCH -n 2
#SBATCH --mem=12G
#SBATCH -t 3-00:00:00
#SBATCH -J vdjtester
#SBATCH -o logs/vdjtester_%j.log
# Launch vdjTester on NAISS Arrhenius.
#   sbatch run_arrhenius.sh [nextflow args]     # full run; the Nextflow head is a small job
#   ./run_arrhenius.sh -stub [args]             # wiring check, directly on the login node
#
# Prerequisites (once):
#   * Cell Ranger unpacked at params.cellranger_dir (conf/arrhenius.config)
#   * nextflow secrets set MI_LICENSE <key>   (stored in ~/.nextflow)
#   * mkdir -p logs   (sbatch writes its log there)
#   * Nextflow at tools/nextflow-${NXF_RELEASE} (the self-contained -dist build):
#       curl -sSfL -o tools/nextflow-26.10.2/nextflow \
#         https://github.com/nextflow-io/nextflow/releases/download/v26.10.2/nextflow-26.10.2-dist
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}"

# The site module is Nextflow 25.10, whose default (v1) config parser lets a
# profile's withName block replace base.config's instead of merging with it.
NXF_RELEASE=26.10.2
module load Java/21.0.11-bdist
export PATH="${PWD}/tools/nextflow-${NXF_RELEASE}:${PATH}"

# Keep Apptainer's own caches in the project, never in $HOME (Apptainer refuses
# images there). The tmpdir needs room for a whole uncompressed image.
export APPTAINER_CACHEDIR="${PWD}/container_cache/apptainer"
export APPTAINER_TMPDIR="${PWD}/container_cache/tmp"
export NXF_APPTAINER_CACHEDIR="${PWD}/container_cache/nextflow"
export NXF_TEMP="${PWD}/container_cache/tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$NXF_APPTAINER_CACHEDIR" logs

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

nextflow run . -profile arrhenius -resume -ansi-log false "${extra[@]}" "$@" 2>&1 | tee "logs/run_$(date +%Y%m%d_%H%M%S).log"
