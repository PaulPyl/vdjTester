Fixtures for stub runs (`-profile test -stub`): empty FASTQs under Cell Ranger's
names, and a samplesheet that mixes local FASTQ directories (`fastq_dir`, relative
to the samplesheet; p1_t1 keeps both libraries in one directory) with ENA runs.

    nextflow run . -profile test -stub --samplesheet assets/test/samplesheet_fastq_dir.csv
