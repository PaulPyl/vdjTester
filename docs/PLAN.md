# vdjTester — consolidation and v0.1 plan

_Written 2026-10-15 after comparing every copy on the laptop and on Dardel._

## 1. What exists today

| copy | git state | last run | notes |
|---|---|---|---|
| **laptop** `~/Projects/Science/vdjTester` | `main`, **2 commits ahead** of GitHub (27 Feb), plus 10 modified files and ~10 untracked files (28 Feb – 2 Mar) | 2 Mar: **compilation error** (`modules/mixcr/main.nf:21`, unescaped `$(pwd)`) | The uncommitted rewrite (opencode, 28 Feb – 2 Mar) does not compile and has structural bugs, listed in §6 |
| **Dardel** `~/Projects/WABI_Common/vdjTester` | `main` at `346cfec` (= GitHub), 1 modified file, 9 untracked | 26 Feb: SRA download OK. 25 Feb: Immcantation, MiXCR and comparison OK | **Has the only end-to-end results.** CellRanger only ran `--version` |
| `~/Projects/Science/vdjTesterOld` | not a git repo | 4 Feb | Earlier prototype. Superseded apart from the `sra-tools_bash` Dockerfile |
| `~/Projects/Science/nextflow_learning/nf-core-vdjtester` | local git, **no remote** | — | nf-core template (5 Mar), with one `downloadcellrangervdjdata` module and an nf-test. Holds the **3.2 GB 10x B-cell example FASTQs** and the VDJ reference |
| GitHub `PaulPyl/vdjTester` | **PUBLIC**, last push 26 Feb | | 4 commits |

**What the 25 Feb Dardel run actually tested:** both MiXCR and Immcantation were fed
CellRanger's *already assembled* `filtered_contig.fasta`, not raw reads. They were then
compared against `BCR_data.tsv`, which comes from a different, older project, not the
same cells. The result: 233 shared sequences, V accuracy 0.987, J 1.0, CDR3 NA. That is
a sanity check of the re-annotation step, not a comparison of the three tools.

## 2. Leaked secret

> **Decision (2026-10-15):** the repo stays public. Paul will apply for a new key when
> the opportunity comes. The current key was **still valid** on 2026-10-15: MiXCR 4.7.0
> runs with it, while a dummy key is rejected by MiLaboratories' licence server. The key
> is now stored with `nextflow secrets set MI_LICENSE` on the laptop and on Dardel, and
> has been removed from every new commit. It is still in the history of commits
> `e2156e9`…`3338b3d`. The MiXCR CLI has no command to list licence usage; ask
> licensing@milaboratories.com or check the account at licensing.milaboratories.com.

- **The MiXCR licence key is committed to `nextflow.config` in all 4 public commits
  (from `e2156e9` onwards).** It is also in the untracked `AGENTS.md`, `nextflow.config`
  on Dardel, and the modified laptop config.
- Treat the key as compromised:
  1. Get a new key from MiLaboratories and invalidate the old one.
  2. Make the repo private now.
  3. Rewrite history with `git filter-repo --replace-text` and force-push. Making the
     repo private does not remove the key from forks or caches; only rotating it does.
- In the future, keep the key out of every file: run `nextflow secrets set MI_LICENSE …`
  on each host and add `secret 'MI_LICENSE'` to the MiXCR process.
- `reference_data/` (`BCR.data_08112023.rds`, `GEX.data_08112023.rds`,
  `BCR_data.tsv`) is gitignored, which is right. Find out which project it came from,
  and keep it out of the repo for good.

## 3. What to keep, and where it goes

| item | where it is now | action |
|---|---|---|
| 6 commits of history | laptop | keep: `main` builds on this |
| uncommitted laptop rewrite, including `cellranger_data/` and `test_lima/` modules and the design docs | laptop | redact the key, then commit it to a branch `wip/2026-03-opencode` as a snapshot. **Do not merge.** Take from it: the 10x download URLs, the `multi` config CSV, the stub blocks, and the useful parts of `Pipeline_Architecture.md` |
| `bin/compare_results.R` (25 Feb; differs from the vdjTesterOld version) | Dardel only | copy into the repo. This script produced the only end-to-end comparison |
| `bin/README_SRA_SUBSAMPLE.md` | Dardel only | copy into `docs/` |
| `bin/compare_results.py` | Dardel = vdjTesterOld (identical) | copy one into the repo |
| `modules/sra/main.nf` (publishDir → storeDir) | Dardel, uncommitted | already superseded by laptop commit `94b1ac0`. Discard |
| `results/` from 25 Feb (about 2 MB) | Dardel | archive to `archive/dardel_2026-02-25/` (gitignored) as a baseline |
| `containers/sra-tools_bash/Dockerfile` | vdjTesterOld | copy to `containers/`. It builds `paulpyl/sra-tools_bash`, which the config references |
| `downloadcellrangervdjdata` module + nf-test, `bin/download_example_data.sh` | nf-core-vdjtester | copy the module in as a reference. Revisit the nf-core layout after v0.1 |
| 10x example FASTQs (3.2 GB) + VDJ reference (12 MB) | nf-core-vdjtester | move to `vdjTester/data/` (gitignored). On Dardel, let the pipeline download its own copy |
| **drop:** CellRanger tarball and unpacked copy (3 GB), `language-server-all.jar`, `bin/nextflow` launcher, `.continue/`, `.vscode/`, `.DS_Store`, `work/`, `.nextflow*`, root Dockerfile/compose | | CellRanger lives outside git on each host. Nextflow comes from pixi (laptop, 26.04.6) or `module load nextflow/26.04.3` (Dardel) |

Delete `vdjTesterOld`, `nf-core-vdjtester` and the old Dardel directory only after the
consolidated repo works on both hosts.

### Done (2026-10-15)

- Uncommitted rewrite → branch `wip/2026-03-opencode` (`651b03d`), key redacted.
- Dardel-only `bin/compare_results.R` and `bin/compare_results.py` → `bin/`. SRA README → `docs/`.
- `containers/sra-tools_bash/Dockerfile` brought over from vdjTesterOld.
- nf-core download module + script → `docs/reference/nf-core-vdjtester/`.
- Dardel 25 Feb results, code and log → `archive/dardel_2026-02-25/` (gitignored).
- 10x FASTQs + VDJ reference moved to `data/`, with symlinks left in nf-core-vdjtester.
- `bin/nextflow` removed; `.gitignore` extended; key removed from `nextflow.config`.
- **Not done yet:** pushing to GitHub; replacing the Dardel checkout; deleting old copies.

## 4. Layout and sync

```
vdjTester/                     git, GitHub (private)
├── main.nf
├── nextflow.config            params + profiles only, no secrets
├── conf/  local.config  dardel.config  test.config
├── modules/local/<tool>/main.nf
├── bin/                       scripts called by processes
├── assets/                    multi-config template, marker gene lists
├── containers/
├── tests/data/                tiny inputs only
├── docs/
└── data/ archive/ results/ work/   gitignored
```

- **Moving code between hosts:** push from the laptop with `git push`, update Dardel
  with `git pull`. Never rsync code.
- **Large files never sync.** Each host fetches its own FASTQs and references through
  the pipeline's download step (with `storeDir`), and has its own CellRanger install.
- **On Dardel:** clone into `~/Projects/WABI_Common/vdjTester`. Put `workDir`,
  `storeDir`, the CellRanger install and the references on project storage, not `$HOME`.
  *Open question: which allocation should this use?*
- **Profiles:**
  - `-profile test,local`: stub runs on the laptop. The Mac is arm64, so real
    containers would run under amd64 emulation; use the laptop for development only.
  - `-profile dardel`: Slurm with `--account`, Apptainer, and escalating `time` on
    retry.

## 5. v0.1: the most basic run that means something

**Goal:** one 10x sample, three tools all starting from the same raw reads, one table
showing per-cell agreement, plus a check of cell identity from gene expression. It
should run end to end on Dardel.

**Data:** 10x "1k B cells, 5′ v2, GEX + BCR" (`sc5p_v2_hs_B_1k_multi_5gex_b`, already
on the laptop).

| step | process | input → output |
|---|---|---|
| 1 | `FETCH_10X_EXAMPLE` (storeDir) | URLs → FASTQs, GEX ref 2020-A, VDJ ref |
| 2 | `CELLRANGER_MULTI` | FASTQs + refs + `assets/multi_config.csv` → filtered feature matrix, `vdj_b/filtered_contig.{fasta,csv}`, `airr_rearrangement.tsv` |
| 3 | `IMMCANTATION` | CellRanger contigs → `AssignGenes igblast` → `MakeDb igblast --10x` → AIRR TSV. Re-annotates CellRanger's assembly; reuse the module that ran on 25 Feb |
| 4 | `MIXCR` | **raw BCR FASTQs** → `mixcr analyze 10x-sc-xcr-vdj` → `exportAirr`. Assembles independently of CellRanger |
| 5 | `COMPARE_VDJ` | three AIRR TSVs → per barcode and locus: chain called (y/n), V/J gene, junction_aa agreement; UpSet-style counts |
| 6 | `EXPRESSION_CHECK` | filtered matrix + all VDJ calls → cell-type × VDJ table (§5a) |
| 7 | `REPORT` | Quarto HTML collecting steps 5 and 6 |

Every tool exports AIRR format: CellRanger writes `airr_rearrangement.tsv`, Immcantation
works in AIRR natively, and MiXCR has `exportAirr`. Standardising on AIRR removes the
need for per-tool column renaming (like the current awk hack in the MiXCR module).
Agreement is measured between tools. The old `BCR_data.tsv` "reference" is dropped,
because it comes from different cells.

### 5a. Expression vs. VDJ

1. **QC and doublets:** filter cells with the usual QC thresholds, then call doublets
   with scDblFinder (R) or scrublet (Python). Record a doublet score per cell; don't
   just drop the doublets.
2. **Marker scores** (module score / AUCell on log-normalised counts):
   - B: `CD79A CD79B MS4A1 CD19 CD22 PAX5 BANK1`
   - T: `CD3D CD3E CD3G CD2 CD247 LCK`
   - Leave out the IG*/TR* V and constant genes, so the expression evidence is
     independent of the VDJ transcripts. Report IG/TR constant-gene expression
     separately as a third, VDJ-adjacent signal.
3. **Classify each cell** as B-like, T-like, both, or neither. Set the score
   thresholds from the data (mixture model or the valley between modes); don't fix
   them in advance.
   - "Both" splits into doublet-flagged cells and **candidate hybrid (dual-expresser)
     cells**, which are "both" but not flagged as doublets.
   - Ambient RNA can also create false "both" calls. Consider SoupX or CellBender, or
     at least report how much ambient signal each marker has.
4. **VDJ class per cell and per tool:**
   - BCR only (IGH and/or IGK/IGL)
   - TCR only (TRA/TRB)
   - both
   - none
   - also whether each chain is productive and paired
5. **Outputs:**
   - Counts table: expression class (B / T / both-doublet / both-hybrid / neither) ×
     VDJ class, per tool.
   - Headline numbers: how many cells look like B or T cells by expression, and what
     fraction of each class gets a VDJ sequence from each tool.
   - Discordant cells listed for inspection: T-like cells with BCR, B-like cells with
     no VDJ.

**This dataset can't answer the T-cell or hybrid questions.** It's a B-cell-enriched
sample with no TCR library, so in v0.1 the T arm only measures B-cell purity. The
second dataset for v0.2 should be a 10x **PBMC 5′ multi with both VDJ-B and VDJ-T
libraries**, ideally with feature barcoding (CD3/CD19 antibodies as protein-level
ground truth for the hybrids). The exact dataset name is still to be chosen from the
10x datasets page. CellRanger multi then needs both `VDJ-B` and `VDJ-T` libraries in
the config, and MiXCR needs one run per chain type, or its combined preset.

### Order of work

1. Fix the secret, make the WIP snapshot branch, collect the files from §3. *(½ day)*
2. Rebuild `main.nf` from the last working commit `346cfec`: one process per tool and
   a named `workflow {}` block per test entry. Get `nextflow run -stub -profile test`
   working locally.
3. On Dardel: `FETCH` + `CELLRANGER_MULTI`. **This is the biggest unknown:** CellRanger
   runs in an `ubuntu` container with a bind mount, and 1k cells needs about 16 cores,
   64 GB, and 1–2 h.
4. `IMMCANTATION` and `MIXCR` on the real outputs. `MIXCR` needs the new licence as a
   Nextflow secret.
5. `COMPARE_VDJ` on AIRR inputs.
6. `EXPRESSION_CHECK` + `REPORT`.
7. v0.2: the PBMC B+T dataset, then subsampling to test how the tools behave at lower
   read depth.

## 6. Known bugs in the uncommitted rewrite (why it's a snapshot, not a base)

- `main.nf` uses names it never includes (`DOWNLOAD_CELLRANGER_DATA`, `RUN_CELLRANGER`,
  `RUN_MIXCR`, `RUN_IMMCANTATION`, `RUN_COMPARISON`), and it has no default unnamed
  workflow.
- Sub-workflows and processes are called several times, once per output:
  `mixcr_workflow(x).ch_tsv`, then `mixcr_workflow(x).ch_per_cell`, and so on. Call
  each once and read `.out.<name>`.
- The `workflow` blocks use `input:`/`output:` and `channel X` instead of
  `take:`/`main:`/`emit:`. `run_cellranger_data` never calls its process.
- Inside the `"""` script blocks: Bash `$(…)` and `$VAR` are not escaped (this is the
  2 March compilation error), and `local` is used outside a function.
- Processes declare both `storeDir` and `publishDir`. The CellRanger stub uses Bash
  `${x:-default}` syntax inside Groovy interpolation.
- The `multi` config CSV hard-codes absolute paths (`/references/…`, `/fastq_samples/…`)
  that don't exist in the task directory.
