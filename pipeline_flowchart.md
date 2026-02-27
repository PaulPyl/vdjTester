# VDJ Tool Evaluation Pipeline Flowchart

## Pipeline Overview
This pipeline evaluates three VDJ (Variable, Diversity, Joining) tools: Immcantation, MixCR, and CellRanger.

## Workflow Structure

```
[Start]
   |
   v
[Input Files]
   |
   v
[Download SRA Fastqs] --> [Run CellRanger]
   |                       |
   v                       |
[Immcantation] <-------- |
   |                       |
   v                       |
[Run MixCR] <----- -- ----|
   |                       |
   v                       |
[Compare Results]
   |
   v
[Output]
```

## Detailed Process Flow

### Main Workflow
1. **Input channels** for fasta, annotations, and reference data are created
2. **SRA Download workflow**:
   - Downloads and subsamples SRA fastqs
   - Stores results in storeDir to avoid re-downloading
3. **CellRanger workflow**:
   - Takes SRA fastqs from download step
   - Processes data with CellRanger tool
4. **Immcantation workflow**:
   - Takes CellRanger output and annotations as input
   - Processes data with Immcantation tool
5. **MixCR workflow**:
   - Takes raw SRA fastqs from download step as input
   - Processes data with MixCR tool
6. **Comparison workflow**:
   - Collects results from all pipelines
   - Compares against reference data
   - Produces final output

### Individual Processes
- **RUN_IMMCANTATION**: 
  - Processes FASTA and annotations 
  - Outputs immcantation.tsv and germ-pass files

- **RUN_MIXCR**:
  - Processes FASTA
  - Outputs multiple files (tsv, per_cell, clns, vdjca, report)

- **RUN_CELLRANGER**:
  - Takes SRA fastqs from download step
  - Processes with CellRanger tool

- **COMPARE_RESULTS**:
  - Takes pipeline results and reference data
  - Produces comparison output

## Container Specifications
- Immcantation: immcantation/suite:4.4.0
- MixCR: ghcr.io/milaboratory/mixcr/mixcr:4.7.0
- CellRanger: ubuntu:22.04
- SRA Download: paulpyl/sra-tools_bash:latest
- Comparison: rocker/tidyverse:latest