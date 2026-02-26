#!/bin/sh
# wrapper script invoked by Nextflow's RUN_CELLRANGER process
set -euo pipefail

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <fastq_dir> [sample_id]" >&2
    exit 1
fi
fastq_dir="$1"
sample_id="${2:-$(basename "$fastq_dir")}"   # default to basename if not provided

# version output
if command -v cellranger &> /dev/null; then
    cellranger --version > cr_version.txt
elif [ -x /opt/cellranger-10.0.0/cellranger ]; then
    /opt/cellranger-10.0.0/cellranger --version > cr_version.txt
else
    echo "ERROR: cellranger binary not found in PATH or /opt" >&2
    exit 1
fi

echo "Running cellranger vdj for sample $sample_id using fastqs in $fastq_dir"
cellranger vdj \
    --id="$sample_id" \
    --fastqs="$fastq_dir" \
    --sample="$sample_id" \
    --localmem=8

echo "Cell Ranger output: $sample_id" > cr_output.txt
