# SRA download and subsample helper

This helper will download SRR runs from SRA (via `prefetch` + `fasterq-dump`) and create small subsampled FASTQ directories suitable for quick testing with Cell Ranger.

Prerequisites
- `sra-tools` (prefetch, fasterq-dump)
- `seqtk` (for subsampling)
- `gzip` / `zcat`

Quick start

Run the example to download and subsample the three default SRRs:

```bash
bin/sra_subsample.sh -n 5000 -o test_fastq
```

This will create `test_fastq/<SRR>_sample/` directories containing files named like `<SRR>_S1_L001_R1_001.fastq.gz` and `_R2_...` which are usable as Cell Ranger inputs.

To supply specific SRR accessions on the command line:

```bash
bin/sra_subsample.sh -n 2000 SRR11233622 SRR11233624
```

Notes
- The script attempts to keep read-pairs synchronized by sampling both R1 and R2 with `seqtk` using the same random seed.
- For full-scale runs, omit `-n` or increase it; downloading full samples may be tens of GB.
