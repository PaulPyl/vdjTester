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
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}"

# The Nextflow module does not bring Java with it.
module load Java/21.0.11-bdist Nextflow/25.10.6-eb

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
