# Candidate test datasets with clonal expansion (vdjTester)

Purpose: find public 10x **5'** datasets with **matched GEX + VDJ libraries** and **raw FASTQs that anyone can download**,
where there is a real clonal signal (dominant malignant B clone, expanded T clones, or age-related loss of diversity).
Complements `docs/datasets.md` (healthy PBMC / naive B datasets).

Research date: 2026-10-16. Nothing was downloaded except metadata (GEO SOFT text, ENA portal API file reports, Synapse REST
folder listings). FASTQ URLs were checked with `curl -sI` only (HTTP 200 + `content-length`).

## How things were verified

- "ENA OK" = the ENA file report lists `fastq_ftp` for the runs and a HEAD request to `https://ftp.sra.ebi.ac.uk/vol1/fastq/...` returned 200
  (done for the specific runs named below; other runs of the same project follow the same URL pattern but were not each HEAD-checked).
- Sizes are ENA `fastq_bytes` (compressed `.fastq.gz`, sum of R1+R2 over all runs of a library).
- ENA URL pattern for `SRR36170840`: `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR361/040/SRR36170840/SRR36170840_{1,2}.fastq.gz`
  (11-character accessions: `SRRxxx/0NN/<acc>/` with `NN` = last two digits of the accession; 10-character ones such as `SRR8315745`
  use `SRR831/005/SRR8315745/`, i.e. the last digit zero-padded). The exact `fastq_ftp` links for any run come from
  `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=<SRR|PRJ>&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv`.
- **None of the SRA/ENA datasets are original 10x tarballs.** They are SRA-regenerated `SRRxxxx_1.fastq.gz` / `_2.fastq.gz` (R1/R2 only, no I1/I2).
  They must be renamed to Cell Ranger style (`<fastq_id>_S1_L001_R1_001.fastq.gz`), and the several runs of one library are the
  lanes/splits (treat as `L001`, `L002`, ...; Cell Ranger only needs unique names). None of the datasets below are 10x BAM uploads.
  Read 1 is not always 26/28 bp (see per dataset); Cell Ranger uses the first 26-28 bp.
- Not verified for any dataset: the number of cells and the clonotype fractions. The "expected signal" statements are from the
  biology/papers, not measured here. A cheap check before committing to a download: the authors' filtered contig annotation CSVs where
  GEO provides them (e.g. `GSE311112_RAW.tar`, 835 MB, HEAD 200 on the NCBI FTP).

## Ranking summary

| Rank | Group | Dataset | Accession | Contrast | Libraries | Raw FASTQ | First-test size | Verdict |
|---|---|---|---|---|---|---|---|---|
| **A1** | A B | SMZL vs healthy spleen/blood | GSE227833 / PRJNA947155 | SMZL (n=3, blood + spleen) vs healthy blood + healthy spleen, same lab | GEX + (BCR, TCR; unannotated) | ENA OK, open | 3-10 GB per sample, 49 GB for everything | **Best matched-control small set** |
| **A2** | A B+T | CLL PBMC on ibrutinib | GSE311112 / PRJNA1368690 | 12 CLL patients, up to 3 time points, ibrutinib-sensitive vs resistant. **No healthy control** | GEX + BCR + TCR | ENA OK, open | **10.6 GB** (cll_1 tp1) | **Best single sample for "dominant clone"** (both B and T) |
| A3 | A B | MM / AL amyloidosis / healthy marrow, sorted B-lineage + plasma cells | GSE293872 / PRJNA1246633 | MM n=11, AL n=5, healthy n=3 (one lacks FASTQs) | GEX + BCR only | ENA OK, open, partly missing | 21-30 GB per sample | Good contrast, but no T cells and big |
| A4 | A B | Follicular lymphoma lymph node (CaClust) | PRJNA1060612, PRJNA1060258, PRJNA1062045, PRJNA1061551, PRJNA978638 | FL only, no control | GEX + BCR-enriched | ENA OK, open | ~47 GB per sample | Backup, single patients |
| **B1** | B T | NSCLC circulating T cells under anti-PD-1/CTLA-4 | GSE266035 / PRJNA1105373 | One patient, pre-treatment + 5 on-treatment points (hashed into 3 GEMs) | GEX + HTO + VDJ-T | ENA OK, open | 12-24 GB per pool, 48 GB total | **Best small T-cell set** |
| B2 | B T | BCC TILs pre/post anti-PD-1 | GSE123813 / PRJNA509910 | Tumour, pre vs post PD-1, 11-14 patients | GEX + TCR (separate libs) | ENA OK, open | 3-20 GB per sample | Solid expansions, tumour only, old 5' v1 |
| B3 | B T | ccRCC tumour + matched blood | GSE121638 (GSE121636 + GSE121637) / PRJNA497905 + PRJNA497906 | 3 patients, tumour vs PBMC | GEX + TCR | ENA OK, open | ~63 GB per sample | Best paired tumour/blood, but large |
| B4 | B T | CTCL (T-cell lymphoma), skin/LN/PBMC + healthy T | GSE290557 / PRJNA1228615 | 35 patients + healthy T | GEX + TCR | ENA OK, open | 80-90 GB per GEX+TCR pair | Strongest T clone, big |
| **C1** | C T | Young vs old vs frail vs cord blood PBMC | GSE157007 / PRJNA659762 | young n=3, old n=6, frail n=5, cord n=3 | GEX + VDJ-T + 4-Ab protein | ENA OK, open | ~75-100 GB per donor (GEX deep) | **Only open ageing set with TCR I could find** |
| C2 | C T+B | Terekhova 2023 blood atlas (166 donors, 25-85 y) | Synapse syn49637038 | age 25-85 | GEX + BCR + TCR + CITE | Synapse, free login, pooled+hashed | very large | Best data, hardest access |

Rejected or access-restricted (not usable for a first test) are listed at the end.

---

## Group A. B-cell clonal expansion

### A1. Splenic marginal zone lymphoma vs healthy donors, 5' scRNA + VDJ (GSE227833)

- Study: Cerapio JP et al., *Integrated spatial and multimodal single-cell transcriptomics reveal patient-dependent cell heterogeneity in splenic marginal zone lymphoma*, J Pathol 2024 (PMID 38828498).
- GEO `GSE227833`, BioProject `PRJNA947155`, SRA `SRP428458`. Title: "A single cell RNA-Seq 5' study of splenic marginal zone lymphomas (SMZL) and spleen healthy donors".
- Samples (8 GSMs, 24 runs, 3 runs per GSM, NextSeq 550, original reads 26+91 bp for blood, ~180 bp spot for spleen):
  - healthy donor 1 **blood** (GSM7109544): SRR23931536, SRR23931537, SRR23931538
  - healthy donor 2 **spleen** (GSM7109543): SRR23931539, SRR23931540, SRR23931541
  - SMZL 1 blood (GSM7109540): SRR23931548, 549, 550; SMZL 1 spleen (GSM7109537): SRR23931557, 558, 559
  - SMZL 2 blood (GSM7109541): SRR23931545, 546, 547; SMZL 2 spleen (GSM7109538): SRR23931554, 555, 556
  - SMZL 3 blood (GSM7109542): SRR23931542, 543, 544; SMZL 3 spleen (GSM7109539): SRR23931551, 552, 553
- Chemistry: Chromium Controller, "Single Cell 5' ... Kit (V3)" as written in GEO (probably 5' v1.1 or v2; not resolved). ~10,000 cells loaded. Frozen cells. Cell Ranger mkfastq + count by the authors.
- Libraries: GEO states the series has GEX, TCR and BCR matrices (`Matrix_SMZL_GEX/T/B`). **The three runs of each GSM are the three libraries, but ENA does not label which is which.**
  By read count the pattern is consistent in every sample: largest run = GEX (e.g. 145 M reads), middle (20-43 M) = most likely BCR, smallest (2-4 M) = most likely TCR. This mapping is an inference, not verified.
  Check it against the GEO supplementary matrices or by `cellranger multi` complaining about VDJ reference mismatches on a tiny test.
- Size (ENA OK, HEAD 200 for SRR23931536, SRR23931537): healthy blood 4.7 GB, healthy spleen 3.3 GB, SMZL blood 4.0-8.7 GB, SMZL spleen 6.8-9.8 GB; **all 24 runs = 48.8 GB**.
- Access: open, no registration. Original FASTQ content via SRA (R1/R2 only).
- Fit: matched controls processed by the same lab, tissue-matched (blood vs blood, spleen vs spleen), tiny, GEX + BCR + TCR in the same cells, so the B/T hybrid check still works.
  Weaknesses: n=3 tumours, controls are two different donors (one blood, one spleen), library-to-run mapping unlabeled, SMZL blood has a lower tumour fraction than CLL blood,
  and cell counts are probably modest (a few thousand per sample).
- Expected signal: per the paper tumour cells form patient-specific B clusters. Expect a single dominant IGH/light-chain clonotype in the SMZL spleens (tumour B cells are the bulk of the spleen B compartment),
  a smaller but visible clone in SMZL blood, and a polyclonal repertoire (top clonotype well under a few percent) in the healthy samples.

### A2. CLL PBMC on ibrutinib, 5' GEX + BCR + TCR (GSE311112)

- Study: Thangavadivel S et al., *Single-Cell RNA-seq Analysis Reveals Distinct Tumor and Immunosuppressive T-Cell Phenotypes in Patients with CLL Treated with Ibrutinib*, Clin Cancer Res 2026 (PMID 42018044). Ohio State (Woyach / Byrd / Blachly groups).
- GEO `GSE311112` (105 GSMs), BioProject `PRJNA1368690`, 216 runs, 716.7 GB total, NovaSeq X Plus, public since 2026-01.
- Design: peripheral-blood PBMC (cryopreserved) from 12 CLL patients, "10X Genomics 5' capture with BCR and TCR sequencing". Ibrutinib-sensitive (IBS, `cll_1`-`cll_4`...) vs resistant (IBR, `cll_5`-`cll_12`).
  3 time points per patient (IBR: diagnosis, BTK C481S detection, relapse; IBS: 3 and 5 years after diagnosis; `timepoint 1` = baseline). Each time point has a `gex`, `bcr` and `tcr` GSM, 2-3 runs each (lane splits).
  **No healthy control in this series.** Use the planned 10x healthy PBMC 10k as the polyclonal comparator (different lab/chemistry, so only a coarse contrast).
- Chemistry: 5' (kit version not stated in GEO); the authors processed with Cell Ranger 3.1.0 (GRCh38). Reads are **2x151** (spot length 302 in ENA), so R1 is 151 bp. Cell Ranger uses the first 26-28 bp, but pass `--chemistry` explicitly if autodetect balks.
- Open, no registration, ENA OK (HEAD 200 verified for SRR36170840 `_1`/`_2` 2.0/2.1 GB, SRR36170838 0.26/0.23 GB, SRR36170836 0.33/0.32 GB).
- Smallest sample sets (sum of all runs of the three libraries):

  | sample | GEX runs | BCR runs | TCR runs | total |
  |---|---|---|---|---|
  | **cll_1 IBS tp1** (GSM9317153/154/155) | SRR36170840, SRR36170841 (8.3 GB) | SRR36170838, SRR36170839 (1.0 GB) | SRR36170836, SRR36170837 (1.3 GB) | **10.6 GB** |
  | cll_2 IBS tp1 (GSM9317162/163/164) | SRR36170822, SRR36170823 (10.2 GB) | SRR36170820, 821 (0.8 GB) | SRR36170818, 819 (1.2 GB) | 12.2 GB |
  | cll_3 IBS tp1 (GSM9317171/172/173) | SRR36170771, 772 (12.7 GB) | SRR36170769, 770 (1.0 GB) | SRR36170767, 768 (2.2 GB) | 15.9 GB |
  | cll_10 IBR tp3 (GSM9317237/238/239) | SRR36170732, 733 (10.7 GB) | SRR36170730, 731 (2.5 GB) | SRR36170728, 729 (2.8 GB) | 16.0 GB |

  Other time points range up to ~50 GB (cll_12 tp1 GEX is 34 GB over 3 runs).
- Cell counts: not checked (the GEO RAW tar, `GSE311112_RAW.tar` 835 MB, contains the authors' matrices/contigs and would give them).
- Fit: the cleanest "one dominant B clone" case: CLL PBMC is mostly the leukaemic clone, so nearly every B cell shares one IGH (+ one IGK or IGL) clonotype.
  The same sample also carries a T compartment with chronic clonal expansions (the paper's T-cell theme: exhausted CD8, Tregs). GEX + BCR + TCR from the same cells, so both lineages are present for the expression-vs-VDJ check.
  Weaknesses: no healthy control, unclear chemistry version, R1 151 bp, SRA renaming needed.
- Expected signal: one IGH clonotype (plus one light chain) in the large majority (likely >50%, often >80-90%) of B cells with productive BCR, naive-like polyclonal B minority; top TCR clones in the percent-to-tens-of-percent range of T cells in some patients.

### A3. Multiple myeloma / AL amyloidosis / healthy bone marrow, sorted B-lineage and plasma cells (GSE293872)

- Study: Gonzalez C et al., *Genetic Evidence against Clonotypic B Lymphocytes as a Reservoir for Plasma-cell Cancers*, Blood Cancer Discov 2026 (PMID 41812162).
- GEO `GSE293872` (38 GSMs), BioProject `PRJNA1246633`, 38 runs, 592 GB, NextSeq 500, public 2026-02.
- Samples: bone marrow aspirates, FACS CD19+CD138- B cells plus CD138+ plasma cells, 10x 5' **v1 and v2** (mixed), Cell Ranger 7.0 (`multi` + `aggr`). 19 pairs of one GEX run + one BCR run:
  MM1-MM11, AL1-AL5, HA1-HA3 ("HA" = healthy adult). Original-style reads: GEX 28+90 (spot 116-119), BCR 2x150 (spot 300).
- **Only B-lineage cells were sorted, so there are no T cells and no VDJ-T library**: the B/T hybrid / TCR side of the pipeline has nothing to find here.
- Gaps in ENA (run exists but 0 FASTQ files): AL1 BCR (SRR33395654), HA2 GEX and BCR (SRR33395643, SRR33395642), MM11 GEX (SRR33395635), MM4 BCR (SRR33395628), MM6 BCR (SRR33395624).
  MM3 GEX (SRR33395631) has a single 26-bp-read file only (11.4 GB, spot length 26): unusable as a normal GEX pair.
- Small sets (ENA OK, HEAD 200 for HA1 and MM8 runs):

  | sample | GEX run | BCR run | total |
  |---|---|---|---|
  | **MM8** | SRR33395621 (16.4 GB) | SRR33395620 (4.6 GB) | 21.0 GB |
  | **AL5** | SRR33395647 (8.2 GB) | SRR33395646 (4.9 GB) | 13.1 GB |
  | **HA1** (healthy) | SRR33395645 (17.1 GB) | SRR33395644 (10.1 GB) | 27.2 GB |
  | MM9 | SRR33395619 (20.3 GB) | SRR33395618 (4.3 GB) | 24.6 GB |
  | MM10 | SRR33395637 (19.6 GB) | SRR33395636 (5.2 GB) | 24.8 GB |
  | HA3 (healthy) | SRR33395641 (32.0 GB) | SRR33395640 (11.3 GB) | 43.3 GB |

  A three-sample first test (HA1 + MM8 + AL5) is about 61 GB.
- Fit: genuine disease/healthy contrast processed by one lab, strong plasma-cell clonality, but B lineage only, healthy n=2 usable, some files missing, sample sizes 13-45 GB, and it is the newest of the candidates (2026).
- Expected signal: CD138+ plasma-cell fraction dominated by one IGH + light-chain clonotype (tumour PCs); healthy marrow PCs and B cells polyclonal.

### A4. Follicular lymphoma lymph nodes with BCR enrichment (CaClust) - backup

- BioProjects `PRJNA1060612` (K6B), `PRJNA1060258` (K7B), `PRJNA1062045` (K4B), `PRJNA1061551` (K5B), `PRJNA978638` (S144); GEO series `GSE252416`, `GSE252344`, `GSE252687`, `GSE252642` (CaClust: linking genotype to transcriptional heterogeneity of FL using BCR and exomic variants).
- Per project: "scRNAseq full transcriptome" (2 runs, ~19+18 GB) + "scRNAseq BCR enriched" (2 runs, ~4.8+4.7 GB) + a PacBio run for K6B. HiSeq 4000, spot length 179 (GEX) and 302 (BCR). ENA OK. ~47 GB per sample.
  PRJNA978638 (S144): GEX one run 30.6 GB + BCR 1.6 GB (32 GB).
- Single patients, no control, chemistry not checked (BCR enrichment implies 5' VDJ). Lymph-node tumour B cells are clonal (FL), so the signal is good, but there is no TCR library and no control. Use only if A1-A3 disappoint.

---

## Group B. T-cell clonal expansion

### B1. Circulating NSCLC T cells under immunotherapy: GEX + HTO + TCR (GSE266035)

- Study: PMID 39863301 ("Divergent gene expression and kinetics of circulating tumor-infiltrating lymphocytes with highly homologous T-cell receptor sequences during immunotherapy"). GEO `GSE266035` (9 GSMs), BioProject `PRJNA1105373`, 18 runs, **48.1 GB**, NovaSeq 6000, public 2025-02.
- Design: one NSCLC patient with response then progression on anti-PD-1/CTLA-4; sorted live CD45+CD3+ PBMC at pre-treatment and 4/8/12/24/32 weeks, hashed with TotalSeq-C hashtags, pooled into 3 GEMs.
  Libraries per pool: 5' GEX, "FB" (the hashtag HTO library) and VDJ-T (TCR). **T cells only, no B cells, no CD19/CD3 protein**: HTO demultiplexing needed.
- Chemistry: Chromium Single Cell 5' kit **v2**, Cell Ranger 6.0.0 + Seurat. Spot length 116 = 26 + 90.
- ENA OK (HEAD 200 verified for SRR28822187 and SRR28822181). Each library has two runs (lane splits):

  | pool | GEX runs | FB (HTO) runs | VDJ runs | total |
  |---|---|---|---|---|
  | BWpool1 (pre, 4wk, 8wk, 12wk) | SRR28822187, SRR28822188 (8.9 GB) | SRR28822175, SRR28822176 (1.4 GB) | SRR28822181, SRR28822182 (1.8 GB) | 12.1 GB |
  | BWpool2 (same time points, 2nd GEM) | SRR28822185, SRR28822186 (9.0 GB) | SRR28822173, SRR28822174 (1.3 GB) | SRR28822179, SRR28822180 (1.7 GB) | 12.0 GB |
  | BWpool3 (24wk, 32wk) | SRR28822183, SRR28822184 (17.2 GB) | SRR28822171, SRR28822172 (3.6 GB) | SRR28822177, SRR28822178 (3.0 GB) | 23.8 GB |

- **Read order caveat (inferred from compressed sizes only):** for SRR28822187 `_1` is 3.1 GB and `_2` is 1.4 GB, i.e. `_1` seems to be the 90-bp cDNA read and `_2` the 26-bp barcode read,
  the opposite of the usual 10x R1/R2. If so, rename `_2 -> R1` and `_1 -> R2`. All other datasets here have the normal orientation (R2 larger, or equal 2x150).
- Fit: smallest T-cell set with modern 5' v2 chemistry; the paper's whole point is expanded circulating T-cell clones with near-identical TCRs, tracked over time (pool1/2 hold pre-treatment vs on-treatment).
  Weaknesses: one patient, no tumour, no B cells (so no B/T hybrid or CD19 check), no control.
- Expected signal: oligoclonal circulating CD8 effector-memory clones; top clonotypes carry a visible share of the T cells (percent-to-tens-of-percent) and expand/contract across time points.

### B2. Basal cell carcinoma TILs pre/post anti-PD-1 (Yost et al. 2019)

- Study: Yost KE et al., *Clonal replacement of tumor-specific T cells following PD-1 blockade*, Nat Med 2019 (PMID 31359002). GEO `GSE123813` (86 GSMs, single-cell), BioProject `PRJNA509910`, SRA `SRP173389`, 86 runs, 650 GB.
- Per sample: one GEX run (`..._RNA`, spot 124 = 26 + 98, HiSeq 4000) and one TCR run (`..._TCR`, 2x150) for su001-su014, pre and post anti-PD-1, several sorted fractions (T cells, tumour CD45, CD39 T cells). **Tumour only, no blood, no B-cell library.**
  Chemistry is early 5' (v1 era, 2018); Cell Ranger 10 compatibility with this chemistry was not tested.
- ENA OK (HEAD 200 verified for SRR8315745 and SRR8315774). Small, clearly-labelled sample pairs:

  | sample | GEX run | TCR run | total |
  |---|---|---|---|
  | su004 pre Tcell | SRR8315745 (4.6 GB) | SRR8315774 (1.7 GB) | 6.3 GB |
  | su004 post Tcell | SRR8315747 (2.3 GB) | SRR8315775 (1.0 GB) | 3.3 GB |
  | su002 pre (All) / post (Tcell) | SRR8315738 (1.2 GB) / SRR8315739 (1.4 GB) | SRR8315770 (0.4 GB) / SRR8315771 (0.7 GB) | 1.6 / 2.1 GB |
  | su014 pre / post Tcell | SRR8890798 (13.1 GB) / SRR8890797 (1.7 GB) | SRR8890815 (6.2 GB) / SRR8890814 (1.0 GB) | 19.3 / 2.7 GB |
  | su009 pre / post Tcell | SRR8890783 (20.3 GB) / SRR8890782 (29.5 GB) | SRR8890800 (3.6 GB) / SRR8890799 (4.2 GB) | 23.9 / 33.7 GB |

  The sample table `GSM -> patient/timepoint/fraction` is in GEO; the `su0xx_*_TCR` libraries correspond to the same-named `*_RNA` library.
- Fit: the paper shows clonally expanded CD8+CD39+ T cells after PD-1 blockade (new clonotypes replacing old ones), so it has a genuine expansion signal and cells from two timepoints.
  Weaknesses: old chemistry, tumour only, T-cell-enriched libraries, some samples tiny.
- Expected signal: a few TCR clones holding a large fraction of the CD8 TILs post-treatment (clonality much higher than blood).

### B3. Clear cell renal carcinoma, tumour + matched PBMC (GSE121638)

- Study: Borcherding N et al., *Mapping the immune environment in clear cell renal carcinoma by single-cell genomics*, Commun Biol 2021 (PMID 33504936).
  GEO super-series `GSE121638`; `GSE121636` (5' GEX) + `GSE121637` (VDJ-T); BioProjects `PRJNA497905` (GEX, 6 runs) and `PRJNA497906` (VDJ, 6 runs); HiSeq 4000.
- 3 patients (GU0700, GU0715, GU0744), tumour (`_T`) and peripheral blood (`_P`), 10x 5'. Per sample: GEX ~55 GB (2x150, 340-360 M reads) + VDJ ~8 GB (SRR8088296-301).
  Total ~380 GB; ~63 GB per sample, ~125 GB per patient (tumour + blood). ENA OK (HEAD 200 verified for SRR8088293).
- Best match for "paired tumour vs blood" (shared clones between tumour and blood). T cells only (TCR), no BCR. GEX is deeper than needed, so pipe-truncate it (recipe at the end).

### B4. Cutaneous T-cell lymphoma (GSE290557) - strongest single T clone, largest data

- GEO `GSE290557`, BioProject `PRJNA1228615` (30 runs, 1.23 TB, NovaSeq 6000; public 2026-03): serial skin, PBMC, lymph node samples from CTCL patients (WU... samples), plus healthy-donor T cells; 5' Gene Expression kit v2, Cell Ranger 3.0.
- 15 `RNA` + `TCR` pairs (e.g. WU1126_RNA2 GSM8816515 SRR32499343 73.8 GB + WU1126_TCR1 GSM8816518 SRR32499340 9.2 GB). GEX deep (2x150, ~500 M reads). ENA OK for listing (HEAD not run).
- A malignant T clone dominates (the malignant CD4 T population shares one TCR), a very strong signal, but a TCR-only, T-lineage-malignancy dataset at ~80 GB per sample pair; only worthwhile with read truncation.

### Also checked for B (not recommended)

- CMV-reactive T cells `GSE235604` / `PRJNA986538`: sorted antigen-specific cells from 13 donors pooled, GEX + VDJ + CSP, 127 GB; not a clean clonality contrast.
- TIL sets with the best-known paired tumour/blood designs are **not openly downloadable**: Wu et al. 2020 Nat (GSE139555, raw at **EGA** EGAS00001003993/94), Zheng et al. 2021 Science (GSE156728, raw at China **GSA**), Liu et al. 2022 (GSE179994, BioProject PRJNA746236: no runs in ENA/SRA).

---

## Group C. Ageing

### C1. Young vs old vs frail vs cord-blood PBMC, 5' GEX + TCR + 4-Ab protein (GSE157007)

- Study: Luo OJ et al., *Multidimensional single-cell analysis of human peripheral blood reveals characteristic features of the immune system landscape in aging and frailty*, Nat Aging 2022 (PMID 37117750).
- GEO `GSE157007` (48 GSMs), BioProject `PRJNA659762`, SRA `SRP279088`, 48 runs, **1.52 TB**, NovaSeq 6000. Fresh PBMC, 10,000 cells loaded, "Chromium Single Cell 5' Feature Barcoding" kits, Cell Ranger 3.1.0.
- Libraries per donor: scRNA (GEX), vdj (TCR) and protein (CCR7, CD45RA, CD4, CD8 antibodies), all 2x150, one run each. **TCR only, no BCR library.**
- Age groups from the GEO sample characteristics:
  - **Young (n=3):** F012, F013, F014
  - **Old (n=6):** F020, F021, F023, OH14, OH15, OH17 (the three `OH` donors have no protein run in ENA)
  - **Frail (n=5):** F002, F004, F006, F007, F008
  - **Cord blood (n=3):** F016, F017, F024
- Runs (ENA OK, HEAD 200 for SRR12535219/233 young F012 and SRR16894504/507 old OH14):

  | donor | group | scRNA | vdj | protein | total |
  |---|---|---|---|---|---|
  | F012 | young | SRR12535219 (81.9 GB) | SRR12535233 (5.7 GB) | SRR12535247 (8.1 GB) | 95.7 GB |
  | F013 | young | SRR12535220 (82.1 GB) | SRR12535234 (7.1 GB) | SRR12535248 (8.1 GB) | 97.3 GB |
  | F014 | young | SRR12535221 (85.1 GB) | SRR12535235 (6.5 GB) | SRR12535249 (7.5 GB) | 99.1 GB |
  | OH14 | old | SRR16894504 (69.4 GB) | SRR16894507 (5.2 GB) | n/a | 74.6 GB |
  | OH15 | old | SRR16894505 (74.7 GB) | SRR16894508 (5.0 GB) | n/a | 79.7 GB |
  | OH17 | old | SRR16894506 (77.0 GB) | SRR16894509 (4.9 GB) | n/a | 81.9 GB |
  | F020 / F021 / F023 | old | SRR12535225 / 226 / 227 | SRR12535239 / 240 / 241 | SRR12535253 / 254 / 255 | 94 / 86 / 91 GB |

- **Too big as is** (~590 M read pairs per GEX library). Practical route for a first test: stream-truncate GEX to ~50 M read pairs (see recipe) and take the VDJ libraries whole: young F012 + old OH14 then costs roughly 20-25 GB of transfer and keeps ~10% of the reads (cell calls are unaffected for 10k loaded cells at 50 M pairs).
  Cord blood is a bonus negative control for expansions (essentially no clones).
- Fit: a real age gradient with TCR + protein and the only open one found; the paper reports the known accumulation of expanded GZMK+ CD8 T cells with age.
  Weaknesses: huge GEX, TCR only (no BCR), n=3 young, protein panel has no CD3/CD19 (only CCR7/CD45RA/CD4/CD8), 5' v1.1 chemistry (Cell Ranger 10 compat not tested).
- Expected signal: young and cord-blood repertoires are near-polyclonal (top clonotype well below 1% of T cells), old donors show higher clonality and expanded CD8 effector-memory/GZMK+ clones (top clones at percent level), frail donors similar or stronger.

### C2. Terekhova et al. 2023 (Immunity), blood atlas of 166 healthy donors aged 25-85

- Terekhova M et al., *Single-cell atlas of healthy human blood unveils age-related loss of NKG2C+GZMB-CD8+ memory T cells and accumulation of type 2 memory T cells*, Immunity 2023.
  317 samples, ~2 M cells, 5' GEX + TCR + BCR + CITE-seq (feature barcode "FB" library) + hashtag multiplexing. Public mirrors exist (AWS Open Data `singlecellhumanbloodatlas`, Hugging Face).
- **Raw FASTQs are on Synapse**, project `syn49637038` (Sage Bionetworks): folders `RAW/`, `BCR_Processed/`, `TCR_Processed/`, `GEX_HTO_processed/`, `Demultiplexing/`, `Cytek/`.
  I listed the folder tree through the anonymous Synapse REST API: `RAW/` has one subfolder per pool (AS044, AS049, AS051-AS062, 14 in total) holding original 10x FASTQs named
  like `ALAW-AS044-1-5pGEX-AS044-1-5pGEX_S1_L004_R1_001.fastq.gz`, `...-BCR-...`, `...-TCR-...`, `...-FB-...` (so they are real bcl2fastq files, not SRA conversions).
  Restriction level `OPEN`, but anonymous users **cannot download** (needs a free Synapse login; file sizes are not visible without it; the project was not downloaded).
- Each pool multiplexes many donors with hashtags, so you cannot fetch one donor in isolation: you need the `Demultiplexing/` outputs and a pool-level run (size unknown, certainly hundreds of GB across libraries).
- Verdict: the best-curated ageing dataset with BCR + TCR, but the access step and pooled design make it unsuitable for a quick first test.

### Not usable for group C

- **Mogilenko et al. 2021, Immunity**: the TCR-seq data are from mouse tissues (GEO `GSE145562`, `GSE154935`, `GSE155006`, mouse). The human PBMC scRNA-seq data (11 young / 10 old healthy men) are only offered through the Artyomov lab browser (`artyomovlab.wustl.edu/immune-aging`), no GEO accession, and that page does not list TCR-seq for them. Cannot be used.
- `GSE231906` (thymus + blood of young/elderly/supercentenarians, scRNA + TCR, PMID 41136751): no runs in ENA/SRA for `PRJNA970218`; GEO only provides processed `GSE231906_RAW.tar`.
- `GSE299953` (Omicron breakthrough infection by age, GEX + TCR + BCR): raw data at China **GSA** (HRA008502), not SRA.
- `GSE136184` (CD8 T cells across the lifespan, 16 donors): 10x **3'** v2, no VDJ.
- Supercentenarian CD4 T cells (PNAS 2019): controlled access (JGA JGAD000324).

---

## Rejected / restricted in Group A (not usable)

- **CLL stereotyped subsets #1/#2/#4**, Karolinska (Oder et al., *Blood Neoplasia* 2026, doi 10.1016/j.bneo.2026.100282): 13 CLL patients, 10x 5' GEX + V(D)J, 138 FASTQ files, ~511 GB.
  Deposited at the Swedish National Data Service (SND `2026-29-1`, doi 10.48723/ep36-6r74) with **"Access to data is restricted"**. Excluded.
- **Follicular lymphoma microenvironment `GSE203610`** (PMID 35687817): 20 FL + 3 reactive lymph nodes, GEX + BCR + TCR, but raw data are at **EGA (EGAS00001006052)**; no runs in SRA/ENA. Excluded.
- **MBL/CLL transition `GSE295489`/`GSE295490`/`GSE295491`** (PMID 41987205; PRJNA1254969/962): 5' RNA + V(D)J BCR in ENA, but each run has a single `SRRxxx.fastq.gz` (ENA layout "PAIRED" yet 90-bp reads), so it looks as if the barcode read is missing or interleaved; not checked further. Multiome RNA libraries from the same series are included. Excluded unless the layout turns out to be fine.
- **10x dataset E donor 3 (ALL patient PBMC, `Donor3_final_Multiplex_fastqs.tar`, 84.0 GB)** from `docs/datasets.md`: HEAD 200 on s3, but the GEX is a `bamtofastq` directory and the config is not a plain multi run; too big.
- Waldenström sets `GSE179221` (CITE-seq) and `GSE235723` (CD19+ marrow, 10x) are GEX-only; no VDJ libraries found. Several myeloma studies (e.g. `PRJEB38620`, `PRJNA1129864` = `GSE271107`, `PRJNA1401834`) were seen in ENA but I did not find BCR libraries or a 5' protocol for them. Not pursued.
- No 10x-hosted dataset with a malignant B or expanded-T clone and open FASTQs was found beyond dataset E (the 10x dataset pages cannot be read from curl, so this was a search-engine check only, not exhaustive).

---

## Recommendation

| Group | Pick | Why |
|---|---|---|
| **A (B-cell clone)** | **GSE227833 (SMZL vs healthy spleen/blood)** for a control contrast, **GSE311112 (CLL), sample `cll_1` tp1** for the strongest clone | A1 has same-lab healthy controls, 5' GEX + BCR + TCR and is only 49 GB in total; A2 is a single 10.6 GB sample with BCR + TCR where a single IGH clone dominates |
| **B (T-cell clone)** | **GSE266035 (NSCLC circulating T cells, pool1, 12 GB)**; backup GSE123813 su004 pre/post (6.3 + 3.3 GB) | smallest modern 5' v2 set with expanding clones over time; BCC is the cleanest tumour-TIL expansion with open FASTQs |
| **C (ageing)** | **GSE157007 (young vs old vs frail vs cord)**, F012 (young) + OH14 (old) with GEX truncated | the only open-FASTQ ageing set with 5' GEX + TCR; Terekhova 2023 is better science but gated and pooled |
| **Overall for the first "expanded clones" test** | **GSE311112, `cll_1` IBS timepoint 1 (GSM9317153/154/155, 10.6 GB)** | one dominant IGH clone in B cells plus a chronically expanded T compartment, GEX + BCR + TCR from the same cells, open, small |

Suggested sequence: (1) `cll_1` tp1 as the "expanded" positive case, against the already-planned 10x healthy PBMC 10k (naive polyclonal comparator);
(2) GSE227833 for same-lab tumour-vs-healthy contrast (HD1 blood vs SMZL1 blood, HD2 spleen vs SMZL1 spleen);
(3) GSE266035 pool1 for T-cell expansion; (4) GSE157007 young vs old for the age contrast.

### Expected signal

- **CLL (A2):** a single IGH heavy-chain clonotype (one light chain, IGK or IGL) in the majority of B cells that have a BCR call, likely >50% and in many patients 80-95%; a minority polyclonal naive-like B fraction.
  The tools should agree on the dominant clonotype, so the metric to compare is the share of the top clonotype among B cells (Cell Ranger vs Immcantation vs MiXCR) and how many cells each tool assigns to it.
  T side: some chronic CD8 clones (top TCR clone a few % to tens of % of T cells).
- **SMZL (A1):** dominant clone in tumour spleens (the paper's patient-specific B clusters), weaker in blood; healthy spleen/blood polyclonal (top clonotype in the low percents or below).
- **MM/AL (A3):** one plasma-cell IGH clone (+ light chain) across most CD138+ cells; healthy marrow polyclonal.
- **NSCLC T cells (B1):** several expanded TCR clones, some shifting across time points; **BCC (B2):** post-PD-1 expanded CD8 clones.
- **Ageing (C1):** clonality increases from cord blood (none) through young to old; old and frail carry more and larger CD8 clones (GZMK+ effector-memory), shown as higher top-clone fraction and lower TCR diversity.
- Numbers above are qualitative expectations; none was measured in this research.

## Practical notes for the pipeline

- **Renaming:** `SRRxxx_1.fastq.gz` / `_2.fastq.gz` -> `<fastq_id>_S1_L00N_R1_001.fastq.gz` / `R2`. Use one `fastq_id` per library and one lane index per SRA run. Exception to check: GSE266035 (probably `_1` = cDNA read, `_2` = barcode read, see B1).
- **Read lengths:** 2x151 (GSE311112, GSE157007, A4, B3, B4, GEX R1 151 bp) is fine for Cell Ranger multi, which only uses the first 26-28 bp of R1 for barcode + UMI for 5' GEX. For 5' VDJ the full R1 is used; for `--r1-length`/`--r2-length` trimming see Cell Ranger docs if mapping quality looks odd.
- **Cell Ranger 10:** `create-bam,false` in `[gene-expression]`, and pass `--chemistry` if auto-detection fails (unknown 5' kit versions for GSE311112 and GSE227833; the 2018 BCC and 2020 frailty sets are 5' v1/v1.1 and not tested with CR 10).
- **Truncating a large GEX download** (untested here, relies on gzip streams being sequential and R1/R2 having identical read order): `curl -s URL_1 | gzip -dc | head -n 200000000 | gzip -c > L001_R1.fastq.gz` and the same for `_2`, giving the first 50 M read pairs without downloading the whole file (curl dies on SIGPIPE).
  Use that for C1 (GEX 70-85 GB), B3 and B4. Keep the VDJ libraries whole: they are small, and truncation would lose clonotypes.
- **HTO demultiplexing** for GSE266035 (4 TotalSeq-C hashtags per pool) needs the FB library and Seurat/`HTODemux` style calling; Cell Ranger multi alone treats it as one sample. Alternatively treat a pool as one sample and ignore the time-point split for a first test.
- **Licences:** GEO/SRA submissions have no explicit open licence; cite the papers. The 10x CC BY 4.0 datasets in `docs/datasets.md` are unaffected.
