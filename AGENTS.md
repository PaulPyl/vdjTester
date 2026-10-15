# VDJ Tool Evaluation Pipeline - Agent Instructions

## Overview
This document provides guidance for agents working on the VDJ Tool Evaluation Pipeline, including recommended tools, documentation URLs, and best practices for Nextflow pipeline development.

---

## Documentation Precedence

### Principle: Documentation as the Source of Truth

When developing or maintaining the VDJ Tool Evaluation Pipeline, always prioritize the [official Nextflow documentation](https://www.nextflow.io/docs/latest/) as the primary reference for:

1. **Syntax and Features**: All Nextflow directives, options, and syntax patterns must align with the latest documentation (version v25.10.4 or newer).
2. **Examples and Code Samples**: Code examples should reflect the actual documentation syntax, particularly:
   - Directives such as `storeDir`, `publishDir`, `cpus`, `memory`, etc., which do not use parentheses per the official docs.
   - Section ordering guidelines and best practices.
3. **Version Compatibility**: Ensure that any examples or configurations leverage features available in the specified Nextflow version.

### Verification Protocol

Before adopting any code example or styling guideline:
- **Reference the official documentation** for the specific feature or directive.
- **Cross-check examples** against the documentation to identify discrepancies in syntax, options, or best practices.
- **Update conflicting examples** to maintain consistency and accuracy.

### Directive Syntax Standards

**storeDir Directive:**
- **Recommended Syntax**: `storeDir "${params.outdir}/salmon"`
- **Avoid**: `storeDir("${params.outdir}/salmon")` (unnecessary parentheses)
- **Purpose**: Directs all generated outputs to the specified directory for efficient storage and cache optimization.

**publishDir Directive:**
- **Recommended Syntax**: `publishDir path: "${params.outdir}/results/salmon", mode: 'copy', saveAs: { it.baseName }, overwrite: true`
- **Purpose**: Publishes output files to a target folder with options for file linking, copying, or hard-linking.
- **Available Modes**: `link`, `copy`, `hardlink`

**Section Ordering:**
- **Standard Sequence**: `tag → storeDir/publishDir → input → output → script → stub → when`
- This ordering ensures logical flow and optimal configuration handling.

---

## Pipeline Structure

### Project Information
- **Project**: VDJ Tool Evaluation Pipeline
- **Configuration**: `documentation-config.json`
- **Nextflow Version**: latest (v25.10.4)
- **Container**: nextflow/nextflow:latest

### Documentation URLs

#### 1. Nextflow Framework

#### Main Documentation
- **Documentation Hub**: https://www.nextflow.io/docs/latest/
- **Get Started**: https://www.nextflow.io/docs/latest/get-started.html
- **Workflow Guide**: https://www.nextflow.io/docs/latest/workflow.html
- **DSL2 Features**: https://www.nextflow.io/docs/latest/dsl2-features.html
- **Executors**: https://www.nextflow.io/docs/latest/nextflow-executor.html

#### Relevant Sections for Pipeline Development
- **Process Stub**: https://www.nextflow.io/docs/latest/process.html#stub
  - Use `stub:` blocks for test data generation with `-stub` flag
  - Stub blocks can be defined before or after the `script:` section
  - Automatic invocation when `-stub` or `-stub-run` CLI options are used

- **Process Script**: https://www.nextflow.io/docs/latest/process.html#script
  - Script section defines executed Bash commands
  - Supports variable substitution with `${variable}` syntax
  - Bash variables should be escaped with `\${variable}` for Nextflow evaluation

- **Section Ordering**: tag → storeDir/publishDir → input → output → script → stub

#### 2. Immcantation

#### Documentation
- **Main Portal**: https://immcantation.readthedocs.io/
- **API Reference**: https://immcantation.readthedocs.io/en/latest/api.html
- **Tutorials**: https://immcantation.readthedocs.io/en/latest/notebooks/tutorials.html
- **Examples**: https://immcantation.readthedocs.io/en/latest/notebooks/examples.html
- **Install Guide**: https://immcantation.readthedocs.io/en/latest/install.html

#### Key Features
- 10x Genomics pipeline integration
- ChangeO analysis tools
- Immcantation visualization
- Human and mouse repertoire analysis
- Container: immcantation/suite:4.4.0

#### Usage in Pipeline
```nf
include { RUN_IMMCANTATION } from './modules/immcantation/main'
```

#### 3. MiXCR

#### Documentation
- **User Guide**: https://mixr.github.io/docs/v4/users_guide/mixr/
- **API Documentation**: https://mixr.github.io/docs/v4/api/mixr/
- **Tutorials**: https://mixr.github.io/docs/v4/tools/
- **Releases**: https://github.com/milaboratory/mixcr/releases

#### Key Features
- Clonotype analysis
- TCR/BCR sequencing analysis
- Single-cell analysis
- License: REDACTED-use-nextflow-secret-MI_LICENSE
- Container: ghcr.io/milaboratory/mixcr/mixcr:4.7.0

#### Usage in Pipeline
```nf
include { RUN_MIXCR } from './modules/mixcr/main'
```

#### 4. CellRanger

#### Documentation
- **Overview**: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
- **Installation**: https://support.10xgenomics.com/single-cell-gene-expression/software/installation-and-system-requirements
- **Pipeline Guide**: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger-count
- **Changelog**: https://support.10xgenomics.com/software/changelog

#### Key Features
- 10x Genomics data analysis
- Single-cell gene expression profiling
- Cell hashing feature
- Integrated pipeline for raw FASTQ analysis
- Container: ubuntu:latest with /opt/cellranger-10.0.0

#### Usage in Pipeline
```nf
include { RUN_CELLRANGER } from './modules/cellranger/main'
```

---

## Recommended Tools

### 1. Nextflow CLI
- **Location**: `bin/nextflow`
- **Usage**: `./bin/nextflow [command]`

#### Commands
```bash
# Lint and validate pipeline
./bin/nextflow lint modules/sra/main.nf

# Run with stub mode
./bin/nextflow run main.nf -stub

# Validate configuration
./bin/nextflow validate
```

### 2. Nextflow Language Server
- **JAR Location**: `bin/language-server-all.jar`
- **Purpose**: Provides IDE features for Nextflow development
- **Supported Extensions**: `.nf`, `.nextflow`, `nextflow.config`

### 3. Apptainer (formerly Singularity)
- **Cache Path**: `/Users/paul/.apptainer/cache`
- **Supported Platform**: `linux/amd64`
- **Configuration**: `nextflow.config`

---

## Process Structure Guidelines

### Section Ordering
For optimal Nextflow process definitions, follow this section ordering:

```nf
process salmon_index {
    // 1. Directives (tag, storeDir, publishDir)
    tag "Salmon Index"

    // storeDir: Directs outputs to ${params.outdir}/salmon
    storeDir "${params.outdir}/salmon"

    // publishDir: Publishes outputs to target directory
    publishDir path: "${params.outdir}/results/salmon",
               mode: 'copy',
               saveAs: { it.baseName },
               overwrite: true

    // 2. Input section
    input:
        path transcriptome

    // 3. Output section
    output:
        path 'index', emit: transcriptomeIndexed
        path 'index/stats.json', emit: stats

    // 4. Script section
    script:
    """
        salmon index --threads $task.cpus -t ${transcriptome} -i index
        cp index/refseq.bin ${transcriptome}.refseq
        cp index/info.json stats.json
    """

    // 5. Stub section (for -stub flag support)
    stub:
    """
        echo "Running in stub mode - Salmon Index"
        mkdir -p index
        touch index/seq.bin
        touch index/refseq.bin
        touch index/info.json
    """

    // 6. Conditional execution
    when: { transcriptome.exists() }
}
```

### Best Practices

#### 1. Stub Mode Support
- Use `stub:` blocks for test data generation
- Stub blocks are automatically invoked with the `-stub` flag
- Generate dummy outputs for downstream tool verification

#### 2. Script Syntax
- Quote file paths in Bash conditionals: `if [ -f "script.sh" ]; then`
- Escape Nextflow variables in Bash: `\${variable}`
- Use triple-quoted strings for multi-line scripts

#### 3. Parameter Management
- Define shared parameters in `nextflow.config`
- Pass parameters through workflow: `${params.parameter}`
- Use subsample parameters for SRA data: `subsample_reads=1000`, `subsample_seed=42`

---

## Development Workflow

### 1. Pipeline Execution
```bash
# Main pipeline execution
./bin/nextflow run main.nf \
    --params-file params.yml \
    -stub

# Lint check
./bin/nextflow lint modules/ -c
```

### 2. Module Import
```nf
// Import reusable modules
include { RUN_IMMCANTATION } from './modules/immcantation/main'
include { RUN_MIXCR } from './modules/mixcr/main'
include { RUN_CELLRANGER } from './modules/cellranger/main'
include { RUN_SRA_DOWNLOAD } from './modules/sra/main'
```

### 3. Test Data Generation
When the `-stub` flag is enabled:
- SRA module generates test data from `test_data/fastq_samples/SRR11233622`
- Creates dummy outputs for CellRanger and MixCR validation
- Copies reference data files for consistency

---

## Configuration Files

### 1. `documentation-config.json`
Contains project metadata, documentation URLs, and environment settings.

### 2. `nextflow.config`
Defines parameters, container settings, and process configurations.

### 3. `opencode.json`
Configures language server integration and supported file extensions.

---

## References

- **Nextflow Documentation**: https://www.nextflow.io/docs/latest/
- **Immcantation**: https://immcantation.readthedocs.io/
- **MiXCR**: https://mixr.github.io/docs/
- **CellRanger**: https://support.10xgenomics.com/single-cell-gene-expression/software/overview
- **GitHub Repository**: [Link to repository]
- **Docker Hub**: https://hub.docker.com/r/milaboratory/mixcr

---

*Last Updated: 2026-02-28*
