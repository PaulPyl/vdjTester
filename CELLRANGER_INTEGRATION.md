# VDJ Tool Evaluation Pipeline - CellRanger Example Data Integration

## Overview

This pipeline has been updated to use **CellRanger's official example data** for VDJ and Gene Expression analysis. The pipeline now downloads 10x Genomics' "Human B cells from a healthy donor" dataset and processes it through multiple bioinformatics tools for comparative evaluation.

## What Changed

### 🆕 New Module: CellRanger Data Download

A new module has been added to handle downloading and preparing the CellRanger example data:

**Location**: `modules/cellranger_data/main.nf`

**Downloads**:
1. **FASTQ Files** (Two libraries):
   - Gene Expression (5' GEX): `sc5p_v2_hs_B_1k_5gex_fastqs`
   - VDJ/BCR Sequencing: `sc5p_v2_hs_B_1k_b_fastqs`

2. **Gene Expression Reference**:
   - Human reference (GRCh38) - 2020-A
   - Used for 5' gene expression analysis

3. **V(D)J Reference**:
   - 5.0 V(D)J human reference with alternate contigs and Ensembl genes
   - Used for BCR/V(D)J analysis

4. **CellRanger Multi Config CSV**:
   - Configuration file specifying pipeline parameters
   - Placed in `reference_data/cellranger_config.csv`

### 🔧 Updated Modules

#### CellRanger Module (`modules/cellranger/main.nf`)
- **Updated** to accept multiple inputs:
  - FASTQ samples
  - Gene Expression reference
  - V(J)D reference
  - Configuration CSV
- **Outputs**:
  - Gene Expression outputs (`ch_gex_outputs`)
  - VDJ outputs (`ch_vdj_outputs`)
  - Web summaries (`ch_web_summary`)
  - Metrics (`ch_metrics`)

#### Immcantation Module (`modules/immcantation/main.nf`)
- **Converted** from process to workflow architecture
- **Inputs**: VDJ data and annotations
- **Outputs**: TSV results and DB-pass files

#### MixCR Module (`modules/mixcr/main.nf`)
- **Converted** to workflow architecture
- **Inputs**: FASTQ samples
- **Outputs**: Multiple channels (TSV, per-cell, clones, VDJCA, reports)

#### Comparison Module (`modules/utils/comparison.nf`)
- **Converted** to workflow architecture
- **Inputs**: Pipeline TSVs and reference data
- **Outputs**: Summary CSV, plots, HTML report

### 🔌 Decoupled SRA Module

The SRA download module has been **decoupled** from the main pipeline:
- Main workflow now uses CellRanger example data
- SRA module remains available for future integration
- Test workflows preserved for modular testing

## Dataset Information

### Example Data: Human B cells from a healthy donor

| Property | Details |
|----------|---------|
| **Sample** | Human B cells from a healthy donor |
| **Cell Count** | ~1,000 cells |
| **Technology** | 10x Genomics Single Cell 5' |
| **Data Type** | Combined 5' GEX + V(D)J Sequencing |
| **Version** | 6.0.0 |

### Data Sources

**Documentation**: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-multi

**Download Sources**:
- **FASTQ**: `https://cf.10xgenomics.com/samples/cell-vdj/6.0.0/...`
- **GEX Reference**: `https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz`
- **VDJ Reference**: `https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz`

## Running the Pipeline

### Main Workflow

```bash
# Run complete pipeline
nextflow run main.nf -profile docker

# Specify output directory
nextflow run main.nf --outdir ./results

# Use local profile for development
nextflow run main.nf -profile local
```

### Parameter Override

```bash
# Custom CellRanger run ID
nextflow run main.nf --cellranger_run_id MyB_Cell_Ana

# Specify output directory
nextflow run main.nf --outdir ./my_results
```

### Test Workflows

```bash
# Test data download only
nextflow run main.nf -w test_cellranger_data

# Test CellRanger pipeline
nextflow run main.nf -w test_cellranger

# Test Immcantation pipeline
nextflow run main.nf -w test_immcantation

# Test MixCR pipeline
nextflow run main.nf -w test_mixcr

# Test comparison analysis
nextflow run main.nf -w test_comparison

# Stub mode (faster testing)
nextflow run main.nf -w test_cellranger -stub
```

## Output Structure

```
results/
├── cellranger_downloads/     # Downloaded FASTQs and references
│   ├── fastq_samples/
│   │   ├── sc5p_v2_hs_B_1k_5gex_fastqs/
│   │   └── sc5p_v2_hs_B_1k_b_fastqs/
│   ├── references/
│   │   ├── refdata-gex-GRCh38-2020-A/
│   │   └── refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0/
│   └── cellranger_config.csv
├── cellranger_output/        # CellRanger multi results
├── immcantation/             # Immcantation results
├── mixcr/                    # MixCR results
└── comparison/               # Comparison analysis results
    ├── comparison_summary.csv
    ├── *.png
    └── comparison_report.html
```

## Configuration File

### cellranger_config.csv

Located in `reference_data/cellranger_config.csv`, this file configures the CellRanger multi pipeline.

## Parameters

### Core Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input_fasta` | `reference_data/filtered_contig.fasta` | Input FASTA file |
| `input_annotations` | `reference_data/filtered_contig_annotations.csv` | Annotation file |
| `reference_data` | `reference_data/BCR_data.tsv` | Reference dataset |
| `outdir` | `results` | Output directory |
| `cellranger_run_id` | `HumanB_Cell_multi` | CellRanger run identifier |

### License Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `mixcr_license` | `E-KUANTQVF...` | MixCR activation license |

### Reserved Parameters (SRA Module)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `subsample_reads` | `1000` | Number of reads for subsampling |
| `subsample_seed` | `42` | Random seed for reproducibility |

## Benefits of This Update

1. ✅ **Standardized Testing**: Uses official 10x Genomics benchmark dataset
2. ✅ **Reproducible Data**: Pre-defined datasets ensure consistent results
3. ✅ **Multi-Tool Comparison**: Evaluates CellRanger, Immcantation, and MixCR
4. ✅ **Modular Architecture**: Each pipeline component tested independently
5. ✅ **Ready for Production**: Full pipeline with comprehensive output

## Future Integration

The SRA download module remains available for:
- Integration with custom SRA datasets
- Large-scale sample processing
- Cloud-based data retrieval

## Documentation Links

- **CellRanger Tutorial**: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-multi
- **Nextflow Documentation**: https://www.nextflow.io/docs/latest/
- **Immcantation Docs**: https://immcantation.readthedocs.io/
- **MixCR Docs**: https://mixr.github.io/docs/

---

*Last Updated: 2024*  
*Pipeline Version: 2.0*  
*CellRanger Version: 10.0.0*