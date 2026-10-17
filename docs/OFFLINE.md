# Running vdjTester without internet (UPPMAX Bianca and similar)

For sensitive human data on an air-gapped cluster: Slurm + Apptainer/Singularity, no
internet on login or compute nodes, files brought in through a gateway (on Bianca, the
wharf). Everything the pipeline would otherwise download is collected beforehand, on a
machine with internet, into one **offline bundle** with a checksum manifest.

Contents:

1. [Data protection: what runs where](#1-data-protection-what-runs-where)
2. [Make the bundle (online)](#2-make-the-bundle-online)
3. [Transfer and verify](#3-transfer-and-verify)
4. [Cell Ranger and MiXCR](#4-cell-ranger-and-mixcr)
5. [Samplesheet with local FASTQs](#5-samplesheet-with-local-fastqs)
6. [Run](#6-run)
7. [Reference: offline parameters and checks](#7-reference-offline-parameters-and-checks)

## 1. Data protection: what runs where

```
 online machine (no sensitive data)           sensitive-data system (e.g. Bianca)
 ----------------------------------           -----------------------------------
 make_offline_bundle.sh                       FASTQs (arrive by the system's own route)
   pipeline code, Nextflow,          wharf    unpack.sh -> pipeline, images, references
   container images, references    ------->   every read-level step:
   (public data only)                           Cell Ranger, MiXCR, TRUST4, Immcantation
                                               every per-cell and per-clone table
                                   <-------   only aggregate outputs, after approval
```

- **Into the system:** only the bundle, which holds public software and public
  references, and no data. FASTQs never pass through the online machine.
- **Inside the system:** everything that touches reads or individual cells. All four
  callers work from reads or from Cell Ranger's contigs (Immcantation re-annotates
  Cell Ranger's assembled contigs), so they all run inside, as does everything
  downstream of them. The Nextflow `work/` directory holds copies of all of it.
- **Out of the system:** only what the project's approval allows. Which outputs are
  aggregate:

| output (under `--outdir`) | content | class |
|---|---|---|
| `cellranger/`, `mixcr/`, `trust4/`, `immcantation/` | contigs, per-cell annotations, count matrices | **read/cell level: stays** |
| `cellranger/*/per_sample_outs/web_summary.html` | QC plus a table of top clonotypes with CDR3 sequences | **sequence level: stays** |
| `compare/*/*.vdj_calls_long.tsv`, `*.vdj_cell_classes.tsv` | one row per cell barcode: genes, CDR3, UMIs | **cell level: stays** |
| `expression/*/*.expression_cells.tsv` | one row per cell: scores, UMAP | **cell level: stays** |
| `clones/*/*.clone_table.tsv`, `*.top_clones.tsv`, `*.clone_tracking.tsv`, `*.dominant_clone_agreement.tsv` | clonotypes by CDR3 sequence | **sequence level: stays** unless the approval covers receptor sequences |
| `vdjtester_report.html` | all aggregate tables, **plus** the largest lineages per sample labelled by CDR3 and the 6 largest CDR3 clonotypes per patient | **sequence level: stays** unless the approval covers receptor sequences |
| `compare/*/*.vdj_tool_summary.tsv`, `*.vdj_pairwise_agreement.tsv`, `*.vdj_v_disagreements.tsv` | counts per tool / tool pair / V-gene pair | aggregate |
| `expression/*/*.expression_summary.tsv`, `*.expression_vs_vdj.tsv`, `*.marker_scores.png` | counts and score distributions | aggregate |
| `clones/*/*.clonality.tsv`, `*.timepoint_overlap.tsv`, `*.clone_summary.tsv`, `*.clone_pairwise.tsv` | diversity indices, overlap and agreement statistics | aggregate |
| `cellranger/*/per_sample_outs/metrics_summary.csv` | sequencing and cell-calling metrics | aggregate |

B- and T-cell receptor sequences are individual-level data (a dominant CDR3 can
identify a person's clone, e.g. a leukaemia), so the report as it stands is not an
aggregate output. Aggregate tables can still carry small cell counts; check them
against the project's disclosure rules before export.

## 2. Make the bundle (online)

On a machine with internet, Nextflow and git, from the repository. References and
images are copied from an existing `store_dir` and image cache where they are there
(the Mac and Dardel runs have both), so the bundle holds exactly what was used there,
and are downloaded or pulled otherwise (pulling needs `apptainer` or `singularity`).

```bash
bin/make_offline_bundle.sh --dry-run --store data/store --container-cache container_cache/nextflow
```

lists every item with its size and source and writes nothing. Then:

```bash
bin/make_offline_bundle.sh --store data/store --container-cache container_cache/nextflow \
    --cellranger-tar ~/Downloads/cellranger-10.0.0.tar.gz /path/to/vdjtester_bundle
```

| option | default | |
|---|---|---|
| `--commit REF` | `HEAD` | pipeline code at this commit (`git archive`, so committed work only) |
| `--reference-set` | `both` | `default`, `imgt` or `both`; `imgt` adds the IMGT reference and MiXCR's IMGT library |
| `--skip-mixcr` | | leave out MiXCR's image and library |
| `--nextflow-version` | local `nextflow -v` | a self-contained Nextflow (`-dist`, needs Java 17+ on the target) |
| `--plugins` | none | Nextflow plugins; the pipeline uses none |
| `--cellranger-tar` | | your own Cell Ranger download (see [4](#4-cell-ranger-and-mixcr)) |
| `--arch` | `amd64` | architecture of images that have to be pulled |

What goes in (sizes from the Mac's store, 2026-10-17; ~18.5 GB in all):

| bundle path | size | source |
|---|---|---|
| `pipeline/vdjTester-<commit>.tar.gz` | 0.1 MB | `git archive` |
| `nextflow/nextflow-<version>-dist` | 40 MB | GitHub releases |
| `containers/immcantation-suite-4.4.0.img` | 1.3 GB | `immcantation/suite:4.4.0` |
| `containers/community.wave.seqera.io-library-python_scanpy_…-dd4b1b2c537a5c1c.img` | 830 MB | Wave image of `containers/python_env.yml` |
| `containers/ghcr.io-milaboratory-mixcr-mixcr-4.7.0.img` | 463 MB | MiXCR |
| `containers/quay.io-biocontainers-trust4-1.1.10--h5814d7d_0.img` | 40 MB | TRUST4 |
| `store/refdata-gex-GRCh38-2024-A.tar` | 15.7 GB (10.7 GB as 10x's `.tar.gz`) | 10x GEX reference |
| `store/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0.tar` | 12 MB | 10x VDJ reference |
| `store/trust4_ref_v1.1.10.tar` | 0.6 MB | TRUST4's `hg38_bcrtcr.fa`, `human_IMGT+C.fa` |
| `store/imgt.tar` | 2 MB | IMGT reference built by fetch-imgt + mkvdjref (`--reference_set imgt`) |
| `store/mixcr_library_imgt.202312-3.sv8.tar` | 0.6 MB | MiLaboratories' IMGT library for MiXCR |
| `tools/cellranger-<version>.tar.gz` | 1-2 GB | only with `--cellranger-tar` |
| `unpack.sh`, `README.txt`, `bundle_info.tsv`, `MANIFEST.sha256` | | |

The IMGT reference has to come from a store: IMGT serves only its current release, so
a new download would differ from the one the comparisons here used. If the dry run
shows `store/imgt.tar` as MISSING, run the pipeline once with `--reference_set imgt`
online, or give `--store` a store that has it.

Images are stored under the names Nextflow's Apptainer/Singularity cache uses, so
Nextflow finds them in `container_cache` and never tries to pull.

## 3. Transfer and verify

Move the whole bundle directory in by the system's route. On Bianca that is the wharf
(see UPPMAX's documentation for the current host name and paths):

```bash
sftp <user>-<project>@bianca-sftp.uppmax.uu.se
sftp> put -r vdjtester_bundle
```

which arrives in `/proj/<project>/nobackup/wharf/<user>/<user>-<project>/`. Then, on
the system:

```bash
cd /proj/<project>/nobackup/wharf/<user>/<user>-<project>/vdjtester_bundle
sha256sum -c MANIFEST.sha256          # each file: OK
./unpack.sh /proj/<project>/vdjtester
```

`unpack.sh` checks the manifest again, then lays out the destination as the offline
profile expects with `--bundle_dir`:

```
/proj/<project>/vdjtester/
  vdjTester/              the pipeline
  nextflow/nextflow       Nextflow
  containers/*.img        container_cache
  store/<name>/           store_dir: references, unpacked
  tools/cellranger        cellranger_dir (link to the unpacked Cell Ranger, if bundled)
```

## 4. Cell Ranger and MiXCR

**Cell Ranger** cannot be redistributed: 10x's licence requires each user to download
it from 10x and accept the terms. The bundle script never downloads it; add your own
download with `--cellranger-tar`, or bring the tarball in separately and unpack it,
or use a Cell Ranger already installed on the system (on UPPMAX, `module spider
cellranger`), and pass its directory as `--cellranger_dir`. It runs on the compute node, not in a
container, as on Dardel. Cell Ranger ships a telemetry tool; it has no network to
reach there, and can be switched off (`cellranger telemetry --help`).

**MiXCR** checks its licence online at every run, at
`licensing-api.milaboratories.com:443` (on Dardel compute nodes it failed with
`License error: ConnectionError`). Without internet it fails unless MiLaboratories has
given you a licence that works offline. Ask them for one, describing the setup
(air-gapped HPC, sensitive data). Until then run with `--skip_mixcr`; the comparison
then covers Cell Ranger, Immcantation and TRUST4. An offline run with MiXCR stops
before starting unless either flag is given:

- `--skip_mixcr`: no MiXCR
- `--mixcr_offline_licence`: you have an offline licence. The key goes in as a Nextflow
  secret in the `NXF_HOME` used there: `nextflow secrets set MI_LICENSE <key>`

## 5. Samplesheet with local FASTQs

Instead of ENA accessions (`runs`), a row can give `fastq_dir`: a directory with
`<fastq_id>_S*_L00*_R1_001.fastq.gz` and `..._R2_001.fastq.gz` (Cell Ranger's names,
as `mkfastq`/`bcl-convert` write them; index reads are ignored). A relative path is
relative to the samplesheet. Several libraries can share one directory.

```csv
sample_id,patient,timepoint,feature_types,fastq_id,runs,fastq_dir
P01_t1,P01,1,Gene Expression,P01_t1_GEX,,/proj/sens2026xxx/fastq/run1
P01_t1,P01,1,VDJ-B,P01_t1_BCR,,/proj/sens2026xxx/fastq/run1
P01_t1,P01,1,VDJ-T,P01_t1_TCR,,/proj/sens2026xxx/fastq/run1
P01_t2,P01,2,Gene Expression,P01_t2_GEX,,/proj/sens2026xxx/fastq/run2
P01_t2,P01,2,VDJ-B,P01_t2_BCR,,/proj/sens2026xxx/fastq/run2
```

Rows are checked before anything runs: each needs `runs` or `fastq_dir` (not both),
and each `fastq_dir` both R1 and R2 files for its `fastq_id`. Two different
directories of one sample must not share a name (Cell Ranger gets them side by side).
`assets/test/samplesheet_fastq_dir.csv` is a stub-test example.

## 6. Run

On Bianca, with the Nextflow driver in `tmux`/`screen` on the login node or in a long,
one-core Slurm job (whichever UPPMAX currently asks for):

```bash
module load java                      # Java 17 or later for Nextflow (module spider java)
D=/proj/<project>/vdjtester
export NXF_OFFLINE=true NXF_HOME=$D/nextflow/home
cd /proj/<project>/vdjtester_run      # launch dir: .nextflow/, logs, work/
$D/nextflow/nextflow run $D/vdjTester -profile bianca \
    --bundle_dir $D --slurm_account <project> \
    --samplesheet samplesheet.csv --outdir results --skip_mixcr \
    -resume
```

Add `--reference_set imgt` for IMGT germlines for every method (needs `store/imgt`, and
for MiXCR its IMGT library). `NXF_OFFLINE=true` keeps Nextflow itself from trying the
network (plugin and pipeline updates).

Check the wiring on the system first, without running any tool (stubs, local executor,
no containers):

```bash
$D/nextflow/nextflow run $D/vdjTester -profile bianca,test -stub \
    --bundle_dir $D --store_dir $D/store --samplesheet samplesheet.csv --skip_mixcr
```

Other air-gapped Slurm systems: `-profile offline` with `--slurm_queue`,
`--slurm_queue_cellranger`, `--max_cpus`/`--max_memory`/`--max_time` and
`--container_engine` set for the site (or a small site config like `conf/bianca.config`).

Booleans: give them as flags (`--skip_mixcr`). In Nextflow 26.04, `--skip_mixcr false`
arrives as the string `"false"`, which counts as true.

## 7. Reference: offline parameters and checks

| parameter | default (offline) | |
|---|---|---|
| `bundle_dir` | | the unpacked bundle; sets the three below |
| `store_dir` | `<bundle_dir>/store` | references, as the download steps would have left them |
| `container_cache` | `<bundle_dir>/containers` | pre-pulled images |
| `cellranger_dir` | `<bundle_dir>/tools/cellranger` | directory with the `cellranger` executable |
| `gex_reference`, `vdj_reference` | from `store_dir` | 10x references, given directly |
| `trust4_reference` | from `store_dir` | directory with `hg38_bcrtcr.fa`, `human_IMGT+C.fa` |
| `imgt_reference` | from `store_dir` (`imgt/vdj_imgt`) | Cell Ranger IMGT reference (mkvdjref output) |
| `mixcr_imgt_library` | from `store_dir` | MiXCR IMGT library; a local path is used as is |
| `slurm_account` | | `sbatch -A` (Bianca: the `sens…` project) |
| `slurm_queue`, `slurm_queue_cellranger` | Bianca: `core`, `node` | partitions |
| `slurm_options` | | anything else for `sbatch` |
| `slurm_mem_per_cpu` | false | `--mem-per-cpu` instead of `--mem` |
| `max_cpus`, `max_memory`, `max_time` | 16, 120 GB, 240 h | every request is capped to these |
| `container_engine` | `apptainer` (Bianca: `singularity`) | |
| `container_options` | | extra options, e.g. `-B /proj` if a container cannot see a path |
| `mixcr_offline_licence` | false | MiXCR runs offline only with this (or `--skip_mixcr`) |

Before anything runs, an offline run stops with a message naming what to stage if:

- a reference is neither given as a path nor in `store_dir` (no download step runs
  offline: `FETCH_TAR`, `FETCH_FILES`, `FETCH_IMGT`, `FETCH_ENA` are all bypassed);
- a samplesheet row has only ENA `runs` and its FASTQs are not in `store_dir/fastq/`;
- a container image is missing from `container_cache` (Nextflow would try to pull it);
- `cellranger_dir` has no `cellranger` (not checked for stub runs);
- MiXCR is on without `--mixcr_offline_licence`.
