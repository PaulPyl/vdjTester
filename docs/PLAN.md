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
vdjTester/                     git, GitHub (public)
├── main.nf
├── nextflow.config            params + profiles only, no secrets
├── conf/  base.config  test.config  dardel.config
├── modules/local/<tool>.nf
├── bin/                       scripts called by processes
├── assets/markers.tsv         B/T marker genes
├── run_dardel.sh
├── docs/
├── attic/v0/                  the Feb 2026 pipeline, kept for reference
└── data/ archive/ results/ work/   gitignored
```

- **Moving code between hosts:** push from the laptop with `git push`, update Dardel
  with `git pull`. Never rsync code.
- **Large files never sync.** Each host fetches its own FASTQs and references through
  the pipeline's download step (with `storeDir`), and has its own CellRanger install.
- **On Dardel:** clone into `~/Projects/WABI_Common/vdjTester`, which links to
  `/cfs/klemming/projects/supr/uppstore2017171/paulpyl/`. **Decided 2026-10-15:** same
  allocation as G_Nordmark_nextflow, i.e. compute on **naiss2026-4-384** (NBIS
  support) and storage in uppstore2017171. Downloads go to
  `…/paulpyl/data/vdjTester_store` and Cell Ranger to
  `…/paulpyl/tools/cellranger-10.0.0`, both outside the repo.
- **Profiles:**
  - `-profile test -stub`: wiring check on the laptop, without containers.
  - `-profile dardel` (via `run_dardel.sh`): Slurm on `shared`, with Cell Ranger alone
    on `main` (whole node). Apptainer, with `time` escalating on retry.

## 5. v0.1: the most basic run that means something

**Goal:** one 10x sample, three tools all starting from the same raw reads, one table
showing per-cell agreement, plus a check of cell identity from gene expression. It
should run end to end on Dardel.

**Data:** 10x "1k B cells, 5′ v2, GEX + BCR" (`sc5p_v2_hs_B_1k_multi_5gex_b`, already
on the laptop).

**Language: Python.** The decision was made on 2026-10-15. The Immcantation steps
used here (`changeo-10x`: IgBLAST + Change-O) are Python too, so R isn't needed
anywhere. All Python steps share the scanpy 1.12 image that G_Nordmark_nextflow uses.

| step | process (file) | input → output | status |
|---|---|---|---|
| 1 | `FETCH_TAR` (`modules/local/fetch_tar.nf`) | URLs → FASTQs, GEX ref 2024-A, VDJ ref 7.1.0 in `store_dir` | written |
| 2 | `CELLRANGER_MULTI` | FASTQs + refs → filtered matrix, `vdj_b/`, `vdj_t/`. Writes its config from `params.libraries` | written; config generation tested with a fake `cellranger` |
| 3 | `IMMCANTATION` | Cell Ranger contigs → `changeo-10x` → AIRR `db-pass` TSV | written (same command as the 25 Feb run) |
| 4 | `MIXCR` | **raw** VDJ FASTQs → `mixcr analyze 10x-sc-xcr-vdj` → per-cell `exportClones` | written; **export field names need checking on the first real run** |
| 4b | `TRUST4` | **raw** VDJ FASTQs (and optionally the GEX library) → `*_barcode_airr.tsv` | **planned**, see below |
| 5 | `COMPARE_VDJ` (`bin/vdj_compare.py`) | all tools' calls → per tool, cell and locus top chain; summary and pairwise V/J/CDR3 agreement | written; tested on 10x's published contigs |
| 6 | `EXPRESSION_CHECK` (`bin/expression_check.py`) | filtered matrix + VDJ classes → B/T/hybrid/doublet classes × VDJ (§5a) | written; tested on synthetic data with planted hybrids |
| 7 | `REPORT` | HTML from steps 5 and 6 | later |

`vdj_compare.py` reads each tool's native output and reduces it to the same columns
itself. That proved more robust than relying on each tool's AIRR export:
- Barcodes lose their `-1` suffix.
- Genes keep only the first call, without the allele.
- CDR3s are compared without the conserved C and W/F, because tools disagree on
  whether those belong to the CDR3.

Agreement is measured between tools. The old `BCR_data.tsv` "reference" is dropped,
because it comes from different cells.

**TRUST4 (planned, decided 2026-10-15).** It's a second assembler that works from raw
reads, alongside MiXCR, so Cell Ranger's assembly is no longer the only one.
- Container: `quay.io/biocontainers/trust4:1.1.10--h5814d7d_0`.
- References: `hg38_bcrtcr.fa` and `human_IMGT+C.fa` from the TRUST4 repository,
  fetched with `FETCH_TAR`-style storage.
- 10x 5′ run: `run-trust4 -u R2 --barcode R1 --barcodeRange 0 15 + --UMI R1 --UMIRange 16 25 + -f hg38_bcrtcr.fa --ref human_IMGT+C.fa`.
  Its `*_barcode_airr.tsv` can go through the Immcantation (AIRR) reader.
- Extra: TRUST4 can also run on the **GEX** library, recovering receptors without VDJ
  enrichment. That gives a third view for the expression-vs-VDJ question.

### 5a. Expression vs. VDJ

1. **QC and doublets:** flag (don't drop) cells failing QC (fewer than 200 genes or
   more than 20 % mitochondrial reads), then call doublets with Scrublet
   (`sc.pp.scrublet`). Record a doublet score per cell; don't just drop the doublets.
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
ground truth for the hybrids).

**Chosen (see `docs/datasets.md`):** 10x "Human PBMC from a Healthy Donor, 10k cells
(v2)", 5′ v2, Cell Ranger 5.0.0.
- Libraries: GEX + antibody capture (CD3, CD19, no CD20) + VDJ-B + VDJ-T.
- About 15.9 GB of FASTQs, and 10x ships a `multi` config.
- 95 barcodes already carry both IG and TR contigs in 10x's own output.

Backup: "5k Human PBMCs, TotalSeq-C TBNK, NextGEM 5′" (Cell Ranger 8.0.0, 33 GB).

For v0.2:
- `params.libraries` gets a `VDJ-T` row and an `Antibody Capture` row, plus the
  `[feature]` reference section.
- MiXCR then runs once per VDJ library; the pipeline already does this.

### Order of work

1. ~~Fix the secret, make the WIP snapshot branch, collect the files from §3.~~ Done 2026-10-15.
2. ~~Rebuild `main.nf`; `nextflow run . -profile test -stub` passes.~~ Done 2026-10-15.
3. **Next, on Dardel:**
   1. Replace the old checkout with a fresh clone.
   2. Unpack Cell Ranger 10.0.0 to `…/paulpyl/tools/` (the tarball is already in the
      old checkout's `bin/`).
   3. `./run_dardel.sh -stub`, then a real run.
   4. Watch Cell Ranger on `main` (64 cores, 200 GB).
4. Fix what the first real run reveals, most likely the MiXCR export field names and
   the changeo-10x output file name.
5. Add `TRUST4`.
6. Add `REPORT`.
7. v0.2: PBMC 10k with B+T+ADT. Use the protein data as ground truth for B/T calls
   and hybrids. Then subsample reads to test how the tools behave at lower depth.

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
