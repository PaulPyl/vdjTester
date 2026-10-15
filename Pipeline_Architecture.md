# VDJ Tool Evaluation Pipeline - Architecture Flowchart



## Pipeline Overview



```
┌────────────────────────────────────────────────────────────────────────────────┐
│                   V D J   T O O L   E V A L U A T I O N   P I P E L I N E      │
│                              Entry Point: main.nf                              │
│                           (Nextflow DSL2 Workflow)                             │
└────────────────────────────────────────────────────────────────────────────────┘


                    ┌─────────────────────────────────────────────────────────┐
                    │                      INPUT LAYER                         │
                    └─────────────────────────────────────────────────────────┘
                              │
                   ┌──────────┼──────────┐
                   │          │          │
    ┌──────────────▼┐  ┌────────▼────────┐  ┌──────────────┐
    │  Input Files  │  │  Reference Data │  │  Parameters  │
    └──────────────┘  └─────────────────┘  └──────────────┘
                              │


┌─────────────────────────────────────────────────────────────────────────────────────┐
│                                                                                    │
│                           MAIN WORKFLOW - main.nf                                  │
│                                                                                    │
│  ┌──────────────────────────────────────────────────────────────────────────────┐ │
│  │                                                                            │   │ │
│  │   ┌──────────────────────────────────────────────────────────────────────┐ │   │ │
│  │   │  CHANNEL LAYERS - Data Flow Between Processes                        │ │   │ │
│  │   │                                                                      │ │   │ │
│  │   │   Input Channels:                                                  │ │   │ │
│  │   │   • ch_input_fasta    → FASTQ/FASTA files                          │ │   │ │
│  │   │   • ch_input_annotations  → CSV annotations                        │ │   │ │
│  │   │   • ch_reference      → TSV reference data                         │ │   │ │
│  │   │   • ch_fastqs         → SRA downloads (channel.empty())            │ │   │ │
│  │   │   • ch_cellranger_out → CellRanger results                         │ │   │ │
│  │   │   • ch_pipeline_tsvs  → Consolidated pipeline outputs              │ │   │ │
│  │   │                                                                      │ │   │ │
│  │   │   Key Connections:                                                 │ │   │ │
│  │   │   ┌────────┐      ┌────────┐      ┌────────┐      ┌────────────┐  │   │ │
│  │   │   │SRA_DL  │ ──→ │ CELLRNGR│ ──→  │IMMCANNTN │ ─────▶COMPARE_RSLTS│  │   │ │
│  │   │   └────────┘      └────────┘      └────────┘    └────────────┘  │   │ │
│  │   │       │               │               │                    ▲        │   │ │
│  │   │       │               │               │                    │        │   │ │
│  │   │       ▼               ▼               ▼                    │        │   │ │
│  │   │   ┌────────┐      ┌────────┐      ┌────────┐      ┌────────────┐  │   │ │
│  │   │   │  MIXCR │ ◀─────└────────┤      └────────┤    │TEST_LIMA   │  │   │ │
│  │   │   │(Parallel  │           │            │    │  │(Lima Test  )│  │   │ │
│  │   │   │ Processing)│           │            │    │  └────────────┘  │   │ │
│  │   │   └────────◀──┘           └────────────┘    └─────────────────────│  │   │ │
│  │   │                                                                                │ │   │ │
│  │   └─────────────────────────────────────────────────────────────────────────────┘│ │
│                                                                                    │ │
└────────────────────────────────────────────────────────────────────────────────────┘ │
                                            │                                          │
                                            │                                          │
                    ┌─────────────────────────┼──────────────────────────────────────────┐
                    │                         │                         │                  │
                    ▼                         ▼                         ▼                  │
                    │                         │                         │                  ▼


┌────────━━━━━┐      ┌──────────────┐      ┌─────────────┐      ┌──────────────┐      ┌─────────────┐
│ SRA_DOWNLOAD│      │  CELLRANGER  │      │ IMMCANTATION│      │    MIXCR     │      │ TEST_LIMA   │
│  (Stub✅)   │      │  (Stub✅)     │      │  (Stub✅)    │      │   (Stub✅)   │      │  (Stub✅)   │
└─────────────┘      └──────────────┘      └─────────────┘      └──────────────┘      └─────────────┘
        │                   │                   │                   │                 │
        ▼                   ▼                   ▼                   ▼                 ▼
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                                INPUT SPECIFICATIONS                                        │
├────────────────────────┬────────────────────────┬────────────────────────────────────────┤
│ SRA_DOWNLOAD           │ CELLRANGER               │ IMMCANTATION                    MIXCR │
├────────────────────────┼────────────────────────┼────────────────────────────────────────┤
│ Input:                │ Input:                   │ Input:                          Input:  │
│ • SRA samples         │ • FASTQ sample dirs       │ • FASTA sequences                 │ FASTQ samples                      │
│ • Subsample params    │   (ch_fastq_stub)          │ • Species annotations             │ • Align params                       │
│   - subsample_reads   │ • CellRanger v10.0       │ • 10x Genomics pipeline           │   - mixcr_license                    │
│   - subsample_seed    │ • Reference: hg38        │   changeo-10x                     │   - Human (10x-sc-xcr-vdj)         │
│ Output:               │ Output:                  │ Output:                     Output:│
│ • FASTQ files         │ • Cell counts            │ • immcantation.tsv              │ • align.vdjca (VDJ alignment)      │
│ • fastq_samples/      │ • raw counts matrix      │ • _germ-pass.tsv (germline)     │ • clones.clns (clonotype library)  │
│ Output: fastq_samples/│   • Gene-cell matrices   │   exportClones()                  │ • exportClones (single-cell)       │
│     SRR11233622       │   • UMI counts           │ Output: tsv, db_pass              │ • report.txt                       │
│                       │                       │   - immcantation.tsv                    │ • mixcr.tsv, mixcr_per_cell.tsv    │
└───────────────────────┴────────────────────────┴────────────────────────────────────┘
                                            │
                                            │
                                            ▼

┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        OUTPUT LAYER - RESULTS                                            │
├────────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                                │
│   ┌──────────────────────────────────┐             ┌────────────────────────────────────┐  │
│   │  CELLRANGER OUTPUT               │             │  COMPARE_RESULTS (Final Analysis)  │  │
│   │  Directory: results/cellranger   │             │                                    │  │
│   │                                ┌▲┐             │ ┌────────────────────────────────┐ │  │
│   │  • cellranger_output/          │ │└─────────────┤  comparison_summary.csv          │ │  │
│   │    └─ sample_id/               │ └────────────▶┤  (Feature comparison & statistics) │ │  │
│   │        ├── raw_feature_count/   │               │ ┌────────────────────────────────┐ │  │
│   │        ├── filtered_count/      │               │ │  comparison_report.html          │ │  │
│   │        ├── cloupe.cloupe        │               │ │  (Interactive visualization)     │ │  │
│   │        ├── web_summary/         │               │ └────────────────────────────────┘ │  │
│   │        └── cellranger_version.txt│               │                                │  │
│   └──────────────────────────────────┘               └────────────────────────────────────┘  │
│                                                                                                │
│   ┌──────────────────────────────────┐             ┌────────────────────────────────────┐  │
│   │  IMMCANTATION OUTPUT             │             │  MIXCR OUTPUT                        │  │
│   │  Directory: results/immcantation │             │  Directory: results/mixcr            │  │
│   │                                ┌▲┐             │                                    │  │
│   │  • immancation/                │ │└─────────────┤  • .clns (Clones Library)          │  │
│   │    └─ 10x-Analysis/            │               │  • .vdjca (Alignment Library)      │  │
│   │        └─ BCR-Analysis/        │               │  • .tsv tables (clonotype data)    │  │
│   │              └─ changeo_db/    │               │                                    │  │
│   │                     └─ tsv     │               │ ┌────────────────────────────────┐ │  │
│   └────────────────────────────────┘                │ │  TEST_LIMA                       │  │
│                                                       │  └─results/test_lima/               │  │
│   ┌──────────────────────────────────┐               │  • lima_test.txt (Apptainer test)│  │
│   │  FINAL COMPACTION LAYER          │               │  • Version compatibility check   │  │
│   └──────────────────────────────────┘               └────────────────────────────────────┘  │
└────────────────────────────────────────────────────────────────────────────────────────┘


┌────────────────────────────────────────────────────────────────────────────────────────┐
│                                    PIPELINE EXECUTION FLOW                               │
│                                                                                                │
│  ┌──────────────────────────────────────────────────────────────────────────────────────┐  │
│  │  PHASE 1: Input Processing                                                             │  │
│  │  ┌────────────────────────────────────────────────────────────────────────────────┐   │  │
│  │  │  params.input_fasta        ───────▶ ch_input_fasta                            │   │  │
│  │  │  params.input_annotations  ─────▶ ch_input_annotations                        │   │  │
│  │  │  params.reference_data     ─────▶ ch_reference                                │   │  │
│  │  │                                                                                │   │  │
│  │  │  Key Parameters:                                                               │   │  │
│  │  │  • input_fasta:       ${projectDir}/reference_data/filtered_contig.fasta      │   │  │
│  │  │  • input_annotations: ${projectDir}/reference_data/filtered_contig_annotations.csv│   │  │
│  │  │  • reference_data:    ${projectDir}/reference_data/BCR_data.tsv               │   │  │
│  │  │  • outdir:            ${projectDir}/results                                   │   │  │
│  │  │  • subsample_reads:   1000                                                     │   │  │
│  │  │  • subsample_seed:    42                                                       │   │  │
│  │  └────────────────────────────────────────────────────────────────────────────────┘   │  │
│  └─────────────────────────────────────────────────────────────────────────────────────────┘  │
│                                            │                                                  │
│                                            ▼                                                  │
│  ┌────────────────────────────────────────────────────────────────────────────────────────┐  │
│  │  PHASE 2: Parallel Tool Execution                                                        │  │
│  │                                                                                         │  │
│  │  ┌─────────────────────────────────┐ ┌────────────────────────────────────────────────┐│  │
│  │  │  PRIMARY PATHWAY                │ │  PARALLEL PATHWAY                                ││  │
│  │  │                                 │ │                                                ││  │
│  │  │  ch_fastqs  ─────▶ RUN_CELLRANGER│◀▶├────────────────────────────────────────────┐ ││  │
│  │  │    ↓                              │ │  RUN_MIXCR                                    │ ││  │
│  │  │  ch_cellranger_out ──────────────▶┤ └──▶ sample_path                              │ ││  │
│  │  │                                  │    ↓                                            │ ││  │
│  │  │  ch_input_annotations ─────────────┤ ──────▶ .vdjca → .clns → .tsv              │ ││  │
│  │  │                                  │    ↓                                            │ ││  │
│  │  │                                  │    └──▶ exportClones (single-cell)           │ ││  │
│  │  └──────────────────────────────────┘    └──────────────────────────────────────────┘ ││  │
│  │                                                                                        ││  │
│  └─────────────────────────────────────────────────────────────────────────────────────────││  │
│                                            │                                                  │  │
└────────────────────────────────────────────────────────────────────────────────────────────┴─┘│
                                            │                                                  │
                                            ▼                                                  │
  ┌────────────────────────────────────────────────────────────────────────────────────────────────┐
  │  PHASE 3: Results Consolidation & Comparison                                                      │
  │                                                                                               │
  │  ┌──────────────────────────────────────┐  ┌────────────────────────────────────────────────┐  │
  │  │  DATA CONSOLIDATION                  │  │  COMPARISON ANALYSIS                             │  │
  │  │                                      │  │                                                │  │
  │  │  ch_pipeline_tsvs:                   │  │  ┌────────────────────────────────────────────┐ │  │
  │  │  ┌────────────────────────────────┐   │  │  │  INPUT TO COMPARE_RESULTS                  │ │  │
  │  │  │  .mix(                         │   │  │  └────────────────────────────────────────────┘ │  │
  │  │  │    RUN_IMMCANTATION.out.tsv,   │   │  │                                                │  │
  │  │  │    RUN_MIXCR.out.per_cell      │   │  │  • ch_pipeline_tsvs   → comparison_summary.csv  │  │
  │  │  │  )                             │   │  │  • ch_reference            → comparison_report   │  │
  │  │  └────────────────────────────────┘   │  │  • Quarto rendering → HTML reports              │  │
  │  │      ↓                                │  │                                                │  │
  │  │  ┌────────────────────────────────┐   │  │  OUTPUT ARTIFACTS:                             │  │
  │  │  │  COMPARE_RESULTS (R/Quarto)    │   │  │  ┌────────────────────────────────────────────┐ │  │
  │  │  └────────────────────────────────┘   │  │  │  1. comparison_summary.csv                 │  │  │
  │  │                                        │  │  │     (Feature comparison stats)            │  │  │
  │  │  ┌────────────────────────────────┐   │  │  ├────────────────────────────────────────────┤ │  │  │
  │  │  │  Input Data:                    │   │  │  │ 2. comparison_report.html                 │  │  │
  │  │  │  • ch_pipeline_tsvs             │   │  │  │     (Interactive Quarto visualization)     │  │  │
  │  │  │  • ch_reference                 │   │  │  ├────────────────────────────────────────────┤ │  │  │
  │  │  └────────────────────────────────┘   │  │  │  │ 3. *.png (Visualizations)                │  │  │
  │  │                                        │  │  │     (Feature plots, quality metrics)       │  │  │
  │  └────────────────────────────────────────┘  │  └────────────────────────────────────────────┘ │  │
  └─────────────────────────────────────────────────────────────────────────────────────────────────┘


┌────────────────────────────────────────────────────────────────────────────────────────┐
│                                    CONTAINER STACK                                       │
│                                                                                                │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────────┐  ┌────────────────────┐ │
│  │ immcantation    │  │  ghcr.io        │  │  rocker/tidyverse   │  │  paulpyl/sra-tools │ │
│  │ suite:4.4.0     │  │  milaboratory/  │  │  (RStudio/Tidyverse)│  │  bash:latest       │ │
│  │                 │  │  mixcr:4.7.0    │  │                     │  │                    │ │
│  │ changeo-10x     │  │                 │  │  - R/Tidyverse      │  │  SRA Toolkit       │ │
│  │ - 10x pipeline  │  │  - align        │  │  - compare_results  │  │  - download_fastq  │ │
│  │ - Species info  │  │  - assemble     │  │  - quarto           │  │  - subsample     │ │
│  └─────────────────┘  └─────────────────┘  └─────────────────────┘  └────────────────────┘ │
│                                                                                                │
│  Container Options: --platform linux/amd64                                                      │
|  Apptainer Support: Enabled (cache-path: ~/.apptainer/cache)                                  │
└────────────────────────────────────────────────────────────────────────────────────────┘


┌────────────────────────────────────────────────────────────────────────────────────────┐
│                           STUB MODE - -stub Flag Support                                 │
│                                                                                                │
│  ┌────────────────┐  ┌────────────────┐  ┌────────────────┐  ┌────────────────┐  ┌────────┐ │
│  │ SRA-Download   │  │ CellRanger     │  │ Immcantation   │  │ MixCR          │  │ TestLima│ │
│  │ (Stub✅)       │  │ (Stub✅)        │  │ (Stub✅)        │  │ (Stub✅)         │  │(Stub✅) │ │
│  └────────────────┘  └────────────────┘  └────────────────┘  └────────────────┘  └────────┘ │
│                                                                                          │
│  Stub Capabilities:                                                                      │
│  ┌─────────────────────────────────────────────────────────────────────────────────────┐ │
│  │ • Automatic invocation with -stub or -stub-run flags                               │ │
│  │ • Test data generation from ${projectDir}/test_data/fastq_samples/               │ │
│  │ • Dummy outputs for tool validation                                                │ │
│  │ • Reference data synchronization (*.tsv, *.csv, *.rds, *.fasta)                    │ │
│  │ • Quick pipeline verification without full execution                               │ │
│  └─────────────────────────────────────────────────────────────────────────────────────┘ │
│                                                                                          │
│  Stub Commands:                                                                          │
│  $ ./bin/nextflow run main.nf -stub                                                       │
│  $ ./bin/nextflow run main.nf -stub-run                                                  │
│                                                                                          │
└────────────────────────────────────────────────────────────────────────────────────────┘

┌──────────────────────────────────────────
│                 OUTPUT DIRECTORY                                         STRUCTURE
└──────────────────────────────────────────

results/
├── fastq_samples/
│   └── SRR11233622/
│       ├── Sample_S1_L001_R1_001.fastq.gz
│       └── Sample_S1_L001_R2_001.fastq.gz
├── cellranger/
│   └── cellranger_output/
│       ├── raw_feature_count/
│       ├── filtered_count/
│       ├── web_summary/
│       └── cellranger_version.txt
├── immcantation/
│   └── immancation/10x-Analysis/BCR-Analysis/
│       ├── immcantation.tsv
│       ├── germ-pass.tsv
│       └── changeo_db/
├── mixcr/
│   ├── align.vdjca (alignment library)
│   ├── clones.clns (clonotype library)
│   ├── mixcr.tsv
│   └── mixcr_per_cell.tsv
├── comparison/
│   ├── comparison_summary.csv
│   └── comparison_report.html
└── test_lima/
    └── lima_test.txt

┌────────────────────────────────────────────────────────────────────────────────────────┐
│                           PROCESS SPECIFICATIONS                                         │
├────────────────────────┬────────────────────────┬────────────────────────────────────────┤
│ Process Name           │ Container             │ Key Directives                             │
├────────────────────────┼────────────────────────┼────────────────────────────────────────┤
│ RUN_SRA_DOWNLOAD       │ paulpyl/sra-tools_bash  │ tag, storeDir, when                      │
│ (SRA Download)         │ :latest                 │                                          │
├────────────────────────┼────────────────────────┼────────────────────────────────────────┤
│ RUN_CELLRANGER         │ ubuntu:latest           │ tag, publishDir, stub                    │
│ (10x Genomics)         │                       │                                          │
├────────────────────────┼────────────────────────┼────────────────────────────────────────┤
│ RUN_IMMCANTATION       │ immcantation/suite      │ tag, publishDir                          │
│ (T/B Cell Analysis)    │ :4.4.0                  │                                          │
├────────────────────────┼────────────────────────┼────────────────────────────────────────┤
│ RUN_MIXCR              │ ghcr.io/milaboratory/   │ tag, storeDir, stub                      │
│ (Clonotype Analysis)   │ mixcr:4.7.0             │                                          │
├────────────────────────┼────────────────────────┼────────────────────────────────────────┤
│ COMPARE_RESULTS        │ rocker/tidyverse:latest │ tag, publishDir                          │
│ (R/Quarto Analysis)    │                       │                                          │
├────────────────────────┼────────────────────────┼────────────────────────────────────────┤
│ TEST_LIMA              │ apptainer               │ tag, publishDir, stub                    │
│ (Environment Testing)  │                       │                                          │
└────────────────────────┴────────────────────────┴────────────────────────────────────────┘


┌────────────────────────────────────────────────────────────────────────────────────────┐
│                              ENVIRONMENT CONFIGURATION                                   │
├────────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                                │
│  Platform Specifics: Linux/amd64, Apptainer enabled, Nextflow DSL2                        │
│  Nextflow Version:  Latest (v25.10.4 or newer)                                             │
│                                                                                                │
│  Configuration Files:                                                                         │
│  • main.nf                    → Main workflow definition                                    │
│  • nextflow.config            → Parameters & container settings                             │
│  • documentation-config.json  → Project metadata & documentation URLs                       │
│                                                                                                │
│  Directives Pattern:                                                                        │
│  ────────────────────────────────────────────────────────────────────────────────────────── │
│  tag → storeDir/publishDir → input → output → script → stub → when                          │
│  ────────────────────────────────────────────────────────────────────────────────────────── │
│                                                                                                │
└────────────────────────────────────────────────────────────────────────────────────────┘
