#!/usr/bin/env bash
set -euo pipefail

# sra_subsample.sh
# Download SRR runs (prefetch + fasterq-dump) and subsample a few thousand read-pairs
# Usage: bin/sra_subsample.sh [-n N_READS] [-s SEED] [-o OUTDIR] [SRR ...]
# If no SRR args given, uses bin/srr_list.txt

N_READS=5000
SEED=42
OUTDIR=fastq

while getopts ":n:s:o:" opt; do
  case ${opt} in
    n ) N_READS=$OPTARG ;;
    s ) SEED=$OPTARG ;;
    o ) OUTDIR=$OPTARG ;;
    \? ) echo "Usage: $0 [-n N_READS] [-s SEED] [-o OUTDIR] [SRR ...]"; exit 1;;
  esac
done
shift $((OPTIND -1))

SRRS=()
if [ $# -gt 0 ]; then
  SRRS=("$@")
else
  if [ -f "$(dirname "$0")/srr_list.txt" ]; then
    mapfile -t SRRS < "$(dirname "$0")/srr_list.txt"
  else
    echo "No SRR IDs provided and bin/srr_list.txt not found." >&2
    exit 1
  fi
fi

mkdir -p "$OUTDIR"

for SRR in "${SRRS[@]}"; do
  echo "Processing $SRR"
  # download (prefetch) then dump paired fastqs
  if ! command -v fasterq-dump &> /dev/null; then
    echo "fasterq-dump (sra-tools) not found in PATH. Install sra-tools or run inside container." >&2
    exit 1
  fi

  prefetch -O "$OUTDIR" "$SRR" || true

  echo "Running fasterq-dump for $SRR"
  fasterq-dump --split-files --outdir "$OUTDIR" --threads 4 "$SRR"

  R1="$OUTDIR/${SRR}_1.fastq"
  R2="$OUTDIR/${SRR}_2.fastq"

  if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    # Accept alternative naming (some sra-toolkit versions use SRRxxxxx_1.fastq etc)
    R1_GZ="$OUTDIR/${SRR}_1.fastq.gz"
    R2_GZ="$OUTDIR/${SRR}_2.fastq.gz"
    if [ -f "$R1_GZ" ] && [ -f "$R2_GZ" ]; then
      gzip -d -c "$R1_GZ" > "$R1"
      gzip -d -c "$R2_GZ" > "$R2"
    else
      echo "fastq files for $SRR not found after fasterq-dump" >&2
      exit 1
    fi
  fi

  # Count reads (pairs) in R1
  total_lines=$(wc -l < "$R1")
  total_reads=$(( total_lines / 4 ))
  echo "Total read pairs in $SRR: $total_reads"

  if [ "$total_reads" -le "$N_READS" ]; then
    echo "Less than or equal to requested reads ($N_READS); keeping full files"
    outdir_sample="$OUTDIR/${SRR}_sample"
    mkdir -p "$outdir_sample"
    gzip -c "$R1" > "$outdir_sample/${SRR}_S1_L001_R1_001.fastq.gz"
    gzip -c "$R2" > "$outdir_sample/${SRR}_S1_L001_R2_001.fastq.gz"
  else
    frac=$(awk -v n=$N_READS -v t=$total_reads 'BEGIN{printf "%.8f", n/t}')
    echo "Sampling $N_READS read pairs (fraction=$frac)"
    outdir_sample="$OUTDIR/${SRR}_sample"
    mkdir -p "$outdir_sample"
    # Use seqtk to sample; keep same SEED for both files so pairs remain synchronized
    if ! command -v seqtk &> /dev/null; then
      echo "seqtk not found in PATH. Install seqtk to subsample fastqs." >&2
      exit 1
    fi
    cat "$R1" | seqtk sample -s$SEED - $frac | gzip > "$outdir_sample/${SRR}_S1_L001_R1_001.fastq.gz"
    cat "$R2" | seqtk sample -s$SEED - $frac | gzip > "$outdir_sample/${SRR}_S1_L001_R2_001.fastq.gz"
  fi

  # cleanup raw (optionally keep originals)
  rm -f "$R1" "$R2"
  echo "Wrote sample fastqs to $outdir_sample/"
done

echo "All done. Subsampled FASTQ directories are under $OUTDIR/"
