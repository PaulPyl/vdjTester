#!/bin/sh
# wrapper for SRA download + subsample used by RUN_SRA_DOWNLOAD
set -eu

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <srr_list> <subsample_reads> <subsample_seed>" >&2
    exit 1
fi

SRR_LIST="$1"
N_READS="$2"
SEED="$3"

mkdir -p fastq_samples

if [ ! -f "$SRR_LIST" ]; then
    echo "SRR list not found: $SRR_LIST" >&2
    exit 1
fi

while read -r SRR; do
    [ -z "$SRR" ] && continue
    echo "Downloading $SRR"
    prefetch -O . "$SRR" || true
    fasterq-dump --split-files --threads 4 --outdir . "$SRR"

    R1="${SRR}_1.fastq"
    R2="${SRR}_2.fastq"
    if [ ! -f "$R1" ]; then R1="${SRR}_1.fastq.gz"; fi
    if [ ! -f "$R2" ]; then R2="${SRR}_2.fastq.gz"; fi

    outdir="fastq_samples/${SRR}_sample"
    mkdir -p "$outdir"

    # Count reads and subsample using awk-based paired sampler
    echo "Subsampling to max $N_READS read pairs"
    
    case "$R1" in
        *.gz)
            total_lines=$(zcat "$R1" | wc -l)
            ;;
        *)
            total_lines=$(wc -l < "$R1")
            ;;
    esac
    total_reads=$(( total_lines / 4 ))
    
    if [ "$total_reads" -le "$N_READS" ]; then
        # Fewer reads than target, just gzip and copy
        echo "  Total reads ($total_reads) <= target ($N_READS), keeping all"
        case "$R1" in
            *.gz)
                zcat "$R1" | gzip > "$outdir/${SRR}_S1_L001_R1_001.fastq.gz"
                zcat "$R2" | gzip > "$outdir/${SRR}_S1_L001_R2_001.fastq.gz"
                ;;
            *)
                gzip -c "$R1" > "$outdir/${SRR}_S1_L001_R1_001.fastq.gz"
                gzip -c "$R2" > "$outdir/${SRR}_S1_L001_R2_001.fastq.gz"
                ;;
        esac
    else
        # Subsample: compute every Nth read using deterministic sampling
        echo "  Subsampling from $total_reads to $N_READS read pairs"
        step=$(( total_reads / N_READS ))
        [ "$step" -lt 1 ] && step=1
        
        # Create awk script to sample paired reads deterministically
        # Keep every Nth read starting from a position based on seed
        offset=$(( (SEED % total_reads) + 1 ))
        
        case "$R1" in
            *.gz)
                # For R1: sample the 4-line records
                zcat "$R1" | awk -v step=$step -v offset=$offset '
                    BEGIN { read_num = 0; rec = 0 }
                    NR % 4 == 1 { 
                        read_num++
                        if ((read_num - offset) % step == 0 && read_num > offset) {
                            rec = NR
                        }
                    }
                    (rec > 0) && (NR >= rec && NR < rec + 4) { print }
                    NR == rec + 3 { rec = 0 }
                ' | gzip > "$outdir/${SRR}_S1_L001_R1_001.fastq.gz"
                
                # For R2: use same indices
                zcat "$R2" | awk -v step=$step -v offset=$offset '
                    BEGIN { read_num = 0; rec = 0 }
                    NR % 4 == 1 { 
                        read_num++
                        if ((read_num - offset) % step == 0 && read_num > offset) {
                            rec = NR
                        }
                    }
                    (rec > 0) && (NR >= rec && NR < rec + 4) { print }
                    NR == rec + 3 { rec = 0 }
                ' | gzip > "$outdir/${SRR}_S1_L001_R2_001.fastq.gz"
                ;;
            *)
                awk -v step=$step -v offset=$offset '
                    BEGIN { read_num = 0; rec = 0 }
                    NR % 4 == 1 { 
                        read_num++
                        if ((read_num - offset) % step == 0 && read_num > offset) {
                            rec = NR
                        }
                    }
                    (rec > 0) && (NR >= rec && NR < rec + 4) { print }
                    NR == rec + 3 { rec = 0 }
                ' "$R1" | gzip > "$outdir/${SRR}_S1_L001_R1_001.fastq.gz"
                
                awk -v step=$step -v offset=$offset '
                    BEGIN { read_num = 0; rec = 0 }
                    NR % 4 == 1 { 
                        read_num++
                        if ((read_num - offset) % step == 0 && read_num > offset) {
                            rec = NR
                        }
                    }
                    (rec > 0) && (NR >= rec && NR < rec + 4) { print }
                    NR == rec + 3 { rec = 0 }
                ' "$R2" | gzip > "$outdir/${SRR}_S1_L001_R2_001.fastq.gz"
                ;;
        esac
    fi

    rm -f "${SRR}_1.fastq" "${SRR}_2.fastq" "${SRR}_1.fastq.gz" "${SRR}_2.fastq.gz"
done < "$SRR_LIST"

echo "Wrote subsampled fastq directories to fastq_samples/"
ls -l fastq_samples || true
