# Candidate test datasets for vdjTester v0.2

Purpose: choose a second public test dataset with **GEX + VDJ-B + VDJ-T from the same cells**, ideally with
**antibody capture (CD3 / CD19 / CD20)** as protein ground truth for the planned expression-vs-VDJ check
(B/T marker scoring, doublet calling, B/T hybrid non-doublet cells, cross-tab against BCR/TCR calls by each tool).

v0.1 dataset for reference: `sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex` (GEX + VDJ-B only, no TCR, no protein).

Research date: 2026-10-15. Nothing was downloaded except small metadata files (configs, feature references,
metrics/contig-annotation CSVs, web summaries) and a few 512-byte tar headers via HTTP Range requests to read tar
structure. No FASTQ content was fetched.

## How URLs were verified

- `OK (HEAD 200)` = `curl -sI` returned HTTP 200 and a `content-length`. Sizes below are that `content-length`.
- 10x dataset **pages** (`www.10xgenomics.com/datasets/...`) return HTTP 429 (Vercel bot checkpoint) to curl and a
  registration form ("welcome" gate) when opened in a browser. They were read via a same-origin `fetch()` from a
  browser session (HTTP 200, page HTML parsed). The file URLs on those pages were then HEAD-checked with curl.
  The registration form was not submitted; the `cf.10xgenomics.com` file URLs themselves need no login.
- Some large tarballs 400 on `cf.10xgenomics.com` (apparently a size limit on the CDN) but are 200 on
  `https://s3-us-west-2.amazonaws.com/10x.files/...`. Those are the URLs 10x links on the page.
- `UNVERIFIED` = not checked or not checkable.

All 10x datasets below are licensed **CC BY 4.0** (stated on each dataset page; cite 10x per their citation guidelines).

## Summary table

| # | Dataset | Chemistry | Cell Ranger (by 10x) | Cells (GEX) | VDJ-B / VDJ-T cells | Protein panel (CD3/CD19/CD20) | FASTQ size | multi config |
|---|---|---|---|---|---|---|---|---|
| A | Human PBMC from a Healthy Donor, 10k cells (v2) | 5' v2 NextGEM, Controller, 26x90 | 5.0.0 | 10,548 | 976 / 5,308 | 17 Abs + 3 isotype ctrls; CD3 yes, CD19 yes, CD20 no | ~15.9 GB | yes |
| B | 5k Human PBMCs stained with TotalSeq-C Human TBNK Cocktail, NextGEM 5' | 5' v2 NextGEM, Chromium X, 28x90 | 8.0.0 | 4,661 | 329 / 2,458 | 9 Abs + IgG1 ctrl; CD3 yes, CD19 yes, CD20 no | ~33.1 GB | yes |
| C | 10k Human PBMCs stained with TotalSeq-C Human TBNK Cocktail, GEM-X 5' | GEM-X 5' v3, Chromium X, 28x90 | 8.0.0 | 9,761 | 470 / 4,967 | same TBNK panel as B | ~74.1 GB | yes |
| D | 10k Human PBMCs, 5' v2.0, Chromium X | 5' v2 NextGEM, 26x90 | 6.1.0 | 9,391 | not extracted | **none** (no FB library) | ~42.8 GB | yes |
| E | Demultiplexing 5' Immune Profiling Libraries Pooled with Hashtags (4 donors) | 5' v2 HT | 7.0.0 | ~60k pooled (per donor ~11k) | not extracted | 130 TotalSeq-C Universal Cocktail incl. CD3, CD19, CD20 + 4 hashtags | >250 GB pooled (see below) | yes |
| F | Integrated GEX, TotalSeq-C and BCR, Chromium Connect, 10k PBMCs | 5' v2 NextGEM Connect | 6.1.2 | 11,075 | BCR only, **no TCR** | TBNK | not checked | yes |

Priority 1-3 from the brief: A and B satisfy all three requirements (GEX + VDJ-B + VDJ-T, CD3 + CD19 protein, 1k-10k cells).
C satisfies them too but is large. D has no protein. E has the best ground truth for doublets but is too large to iterate on.
F is unsuitable (no TCR library).

---

## A. Human PBMC from a Healthy Donor, 10k cells (v2)  (10x 5.0.0, Nov 2020)

- Pages:
  - Count/vdj version with all FASTQs: https://www.10xgenomics.com/datasets/human-pbmc-from-a-healthy-donor-10-k-cells-v-2-2-standard-5-0-0  (read via browser fetch, 200)
  - "multi" version (processed outputs + multi config): https://www.10xgenomics.com/datasets/human-pbmc-from-a-healthy-donor-10-k-cells-multi-v-2-2-standard-5-0-0  (read via browser fetch, 200; this page lists **no FASTQs**, use the sibling page's)
- Sample: PBMCs, one healthy male donor aged 27 (AllCells). 16,520 cells loaded, 91% viable, 10,548 cells called in GEX.
- Chemistry: Chromium Next GEM Single Cell 5' Reagent Kits v2 (Dual Index) with Feature Barcode technology (CG000330 Rev A); NovaSeq 6000; paired-end 26x90, i7 10, i5 10 for every library. Chromium Controller.
- Cell Ranger used by 10x: 5.0.0 (`multi` plus separate `count`/`vdj`). GEX `--expect-cells=10000`, ref GRCh38-2020-A, VDJ ref `vdj_GRCh38_alts_ensembl-5.0.0`.
- Libraries: Gene Expression, Antibody Capture ("Feature Barcode"), VDJ-B, VDJ-T. Two lanes each (`lanes 1|2` in the config).
- Cells per library (from the multi web summary / contig annotations): GEX 10,548; VDJ-T 5,308 cells; VDJ-B 976 cells.
  Reads: VDJ-T 35.2 M reads (6.6k reads/cell), VDJ-B mean ~39.8k reads/cell, GEX ~60.5k mean reads/cell.
- Antibody panel (`sc5p_v2_hs_PBMC_10k_feature_ref.csv`, OK 200, 1,550 B), BioLegend TotalSeq-C, all `Antibody Capture`, R2, pattern `^NNNNNNNNNN(BC)NNNNNNNNN`:
  **CD3**, **CD19**, CD45RA, CD4, CD8a, CD14, CD16, CD56, CD25, CD45RO, PD-1, TIGIT, CD127, CD15, CD197 (CCR7), HLA-DR,
  isotype controls IgG1, IgG2a, IgG2b. **No CD20.** CD3 and CD19 are present, which is what the hybrid check needs.
- Quick plausibility numbers (computed from 10x's published `*_vdj_{b,t}_filtered_contig_annotations.csv`):
  976 B barcodes, 5,308 T barcodes, **95 barcodes carry both IG and TR contigs** (about 10% of B-called cells). Those are
  exactly the BCR+TCR hybrid/doublet candidates the planned check is designed to sort out. No TR contigs appear in the B library and no IG contigs in the T library, i.e. no cross-library contamination at contig level.
- FASTQ tarballs (all OK, HEAD 200, `last-modified` 2020-11-16):
  - GEX + FB: https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_f_fastqs.tar  (9,584,926,720 B = 9.6 GB)
  - VDJ-B: https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_b_fastqs.tar  (3,337,297,920 B = 3.3 GB)
  - VDJ-T: https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_t_fastqs.tar  (2,995,712,000 B = 3.0 GB)
  - Total about **15.9 GB** (15,917,936,640 B). Smallest of the candidates with all four libraries.
  - Tar layout, read from tar headers by Range request (not a download): `..._f_fastqs.tar` has top-level directory `sc5p_v2_hs_PBMC_10k_5fb_fastqs/` with files `sc5p_v2_hs_PBMC_10k_5fb_S1_L00{1,2}_{R2,I1,...}_001.fastq.gz`;
    `..._b_fastqs.tar` is `sc5p_v2_hs_PBMC_10k_b_fastqs/sc5p_v2_hs_PBMC_10k_b_S1_L00{1,2}_{R1,R2,I1,I2}_001.fastq.gz`.
    The GEX directory `sc5p_v2_hs_PBMC_10k_5gex_fastqs/` is implied by the 10x config paths but its tar header was not reached (range reads timed out); UNVERIFIED but low risk.
    The files are named with `_S1_L00x_`, so `fastq_id` = `sc5p_v2_hs_PBMC_10k_5fb` etc. matches Cell Ranger's prefix rule.
- 10x multi config (OK 200, 756 B): https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_config.csv
  - libraries: `sc5p_v2_hs_PBMC_10k_5gex` (gene expression), `..._5fb` (antibody capture), `..._b` (vdj-b), `..._t` (vdj-t), `lanes 1|2`.
  - feature ref also at https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_feature_ref.csv (OK 200, identical 1,550 B).
- Other useful files (all OK 200): multi web summary `.../sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_web_summary.html`; per-library `*_vdj_{b,t}_filtered_contig_annotations.csv` (0.42 MB, 2.2 MB); `*_count_filtered_feature_bc_matrix.h5` (29.9 MB) as a cheap reference for expected cell calls; `*_vdj_reference.tar.gz` (42 KB).
- Caveats for Cell Ranger 10.0:
  - The 10x config targets 5.0.0 references. For CR 10.0 use `refdata-gex-GRCh38-2024-A` (https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz, OK 200, 11.5 GB) and
    `refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0` (https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0.tar.gz, OK 200, 3.4 MB). Using 2020-A / vdj 5.0.0 also still works and keeps v0.1 and v0.2 comparable (v0.1 used 2020-A and 5.0.0 VDJ per the 10x tutorial).
  - Since Cell Ranger 8, `create-bam` is mandatory in `[gene-expression]` of the multi config (the 10x 8.0.0 configs below show `create-bam,true`). Add `create-bam,false` to avoid BAMs.
  - Read 1 is 26 bp (16 bp barcode + 10 bp UMI); this is the 5' v2 layout and works with CR 10.0 unchanged.
  - Old FASTQ naming (`S1_L001`) is fine; bcl2fastq-style names.

## B. 5k Human PBMCs stained with TotalSeq-C Human TBNK Cocktail, Chromium NextGEM Single Cell 5'  (10x 8.0.0, Mar 2024)

- Page: https://www.10xgenomics.com/datasets/5k-human-pbmcs-stained-with-totalseq-C-human-TBNK-cocktail-NextGEM  (browser fetch 200)
- Sample: PBMCs, healthy male donor 18-35 (AllCells); 2 M cells stained with BioLegend TotalSeq-C Human TBNK Cocktail (Cat# 399903). Fresh frozen.
- Chemistry: Chromium Single Cell 5' v2 (Dual Index) with Feature Barcoding (CG000330), Chromium X, NovaSeq 6000, about 75,000 read pairs/cell; 28 / 10 / 10 / 90.
- Cell Ranger by 10x: 8.0.0 `multi`, GRCh38-2024-A, `vdj_GRCh38_alts_ensembl-7.1.0`.
- Cells: **4,661** (GEX = antibody = same). VDJ-B 329 cells (225 with productive pair), VDJ-T 2,458 cells (2,117 with productive pair). 68 barcodes carry both IG and TR contigs (about 21% of B-called cells; computed from published filtered contig annotations).
  Mean GEX reads/cell 54,568; median genes/cell 2,296; antibody 1,718 reads/cell, 841 usable.
- Antibody panel (`..._count_feature_reference.csv`, OK 200, 767 B): **CD3, CD19**, CD4, CD8, CD11c, CD14, CD16, CD56, CD45, IgG1 control. **No CD20.**
  The config's feature reference is called `TotalSeqB_TBNK_Panel_Nov2020.csv` (sic, a naming slip in 10x's file; the barcodes are the TotalSeq-C sequences, e.g. CD3 `CTCATTGTAACTCCT`, CD19 `CTGGGCAATTACTCG`, matching dataset A).
- FASTQ tarball (OK, HEAD 200 on both hosts, `last-modified` 2024-03-11):
  https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/5k_Human_PBMC_TotalSeqC_5p_nextgem_Multiplex/5k_Human_PBMC_TotalSeqC_5p_nextgem_Multiplex_fastqs.tar  (33,054,392,320 B = **33.1 GB**)
  (same file at https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/8.0.0/5k_Human_PBMC_TotalSeqC_5p_nextgem_Multiplex/5k_Human_PBMC_TotalSeqC_5p_nextgem_Multiplex_fastqs.tar).
  Per the config the tar unpacks to `5k_Human_PBMC_TotalSeqC_5p_nextgem_fastqs/{gex,ab,vdj_b,vdj_t}`, 4 lanes for GEX/VDJ, two GEX fastq_ids (`_gex1`, `_gex2`). Internal layout from the config only; tar headers not read, UNVERIFIED.
- multi config (OK 200, 1,128 B): https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/5k_Human_PBMC_TotalSeqC_5p_nextgem_Multiplex/5k_Human_PBMC_TotalSeqC_5p_nextgem_Multiplex_config.csv
  Already CR 8+ style (`create-bam,true`), uses GRCh38-2024-A and VDJ 7.1.0, i.e. the closest to a Cell Ranger 10.0 setup without edits.
- Feature reference: https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/5k_Human_PBMC_TotalSeqC_5p_nextgem_Multiplex/5k_Human_PBMC_TotalSeqC_5p_nextgem_Multiplex_count_feature_reference.csv
- VDJ reference tarball (OK 200, 36.6 KB): `..._Multiplex_vdj_reference.tar.gz` in the same directory.
- Downside: 2x the download of A and fewer B cells (329 vs 976). Upside: modern Cell Ranger/reference versions, higher B/T overlap (more hybrid candidates in relative terms).

## C. 10k Human PBMCs stained with TotalSeq-C Human TBNK Cocktail, Chromium GEM-X Single Cell 5'  (10x 8.0.0, Mar 2024)

- Page: https://www.10xgenomics.com/datasets/10k-human-pbmcs-stained-with-totalseq-C-human-TBNK-cocktail-GEM-X  (browser fetch 200)
- Same donor description and TBNK panel as B, but GEM-X 5' v3 (CG000734), Chromium X, about 62,000 read pairs/cell.
- Cells: **9,761**; VDJ-B 470 cells, VDJ-T 4,967 cells; 23 barcodes with both IG and TR contigs. Mean GEX reads/cell 56,536, median genes/cell 3,598.
- FASTQs: https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/8.0.0/10k_Human_PBMC_TotalSeqC_5p_gemx_Multiplex/10k_Human_PBMC_TotalSeqC_5p_gemx_Multiplex_fastqs.tar
  (OK HEAD 200, 74,086,768,640 B = **74.1 GB**). The `cf.10xgenomics.com` form of this URL returns 400; use the s3 host.
- multi config (OK 200): https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/10k_Human_PBMC_TotalSeqC_5p_gemx_Multiplex/10k_Human_PBMC_TotalSeqC_5p_gemx_Multiplex_config.csv
  feature ref (OK 200): `..._count_feature_reference.csv` in the same directory. Same 10-feature panel (CD3, CD19, CD4, CD8, CD11c, CD14, CD16, CD56, CD45, IgG1).
- Too large for iteration; useful later as a "scale-up" or a second chemistry (GEM-X) check.

## D. 10k Human PBMCs, 5' v2.0, Chromium X  (10x 6.1.0, Aug 2021)  - no protein

- Page: https://www.10xgenomics.com/datasets/10-k-human-pbm-cs-5-v-2-0-chromium-x-2-standard-6-1-0  (browser fetch 200)
- GEX + VDJ-T + VDJ-B, 9,391 cells, healthy male donor 30-35, 5' v2 NextGEM, 26x90, ~30k GEX / 15k TCR / 25k BCR reads/cell, `--expect-cells=10000`.
- **No Feature Barcode library**, so no CD3/CD19 ground truth. Not recommended.
- FASTQs: https://cf.10xgenomics.com/samples/cell-vdj/6.1.0/10k_PBMC_5pv2_nextgem_Chromium_X_Multiplex/10k_PBMC_5pv2_nextgem_Chromium_X_Multiplex_fastqs.tar  (OK 200, 42,818,549,760 B = 42.8 GB)
- multi config (OK 200, 1,060 B): `..._Multiplex_config.csv` in the same directory (uses `ref` keys, `lanes any`, `physical_library_id`; needs updating to `reference` for CR 10).
- There is a "without intronic reads" sibling listed by search (`10k-human-pbmcs-5-v2-0-chromium-x-without-intronic-reads-2-standard`), UNVERIFIED, same libraries.

## E. Demultiplexing and Analyzing 5' Immune Profiling Libraries Pooled with Hashtags  (10x 7.0.0, Sep 2022)  - best ground truth, too big

- Page: https://www.10xgenomics.com/datasets/5-hashing-example-with-tabs-2-standard  (browser fetch 200)
- Four donors (2 healthy PBMC, 1 ALL bone marrow, 1 ALL PBMC) each stained with the **TotalSeq-C Human Universal Cocktail (130 Abs; CD3, CD19, **CD20**, CD4, CD8, CD14, CD56 ... checked in the feature reference, 137 data rows)** plus one of four hashtag antibodies, pooled, 5' v2 HT.
  Libraries: GEX (44,856 reads/cell), Antibody Capture (8,052), BCR (5,837), TCR (5,043). Cell Ranger 7.0.0; GRCh38-2020-A, VDJ 7.1.0 (page text) / 7.0.0 (config).
- Hashtags give **cross-donor doublet ground truth** (two hashtags in one droplet), the ideal reference for the B/T doublet question. It is also the only candidate with CD20.
- But: pooled FASTQs are enormous (all HEAD 200 on s3): `5p_hashing_demux_gex_fastq.tar` 184.1 GB, `_ab_fastq.tar` 27.9 GB, `_tcr_fastq.tar` 20.5 GB, `_bcr_fastq.tar` 11.6 GB (total about 244 GB);
  per-donor tars `Donor1_final_Multiplex_fastqs.tar` 91.7 GB, `Donor2_final_Multiplex_fastqs.tar` 79.5 GB. Base: https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/7.0.0/ . Disease donors (ALL) complicate B/T biology (Donor 3 has zero median TRA UMIs).
  multi configs for the per-donor runs exist (e.g. `.../Donor1_final_Multiplex/Donor1_final_Multiplex_config.csv`, OK 200) and use `bamtofastq` GEX; not a plain `cellranger multi` from FASTQs.
- Suggested use: a later (v0.3+) stress test, not v0.2.

## F. Integrated GEX, TotalSeq-C, and BCR Analysis of Chromium Connect Generated Library from 10k Human PBMCs  (10x 6.1.2)  - rejected

- Page: https://www.10xgenomics.com/datasets/integrated-gex-totalseq-c-and-bcr-analysis-of-chromium-connect-generated-library-from-10k-human-pbmcs-2-standard  (browser fetch 200)
- GEX + BCR + Antibody Capture (TBNK), 11,075 cells, but **no TCR library** (same limitation as v0.1). Rejected.

## Older / other 10x sets considered

- PBMCs of a Healthy Donor - 5' Gene Expression with a Panel of TotalSeq-C Antibodies (10x 3.0.0, 2018, `vdj_v1_hs_pbmc2_5gex_protein`; 8,258 cells; CD3, CD19 + 15 others):
  https://www.10xgenomics.com/datasets/pbm-cs-of-a-healthy-donor-5-gene-expression-and-cell-surface-protein-1-standard-3-0-0 (browser fetch 200).
  FASTQ https://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_hs_pbmc2_5gex_protein/vdj_v1_hs_pbmc2_5gex_protein_fastqs.tar (OK 200, 42.9 GB), feature ref OK 200.
  Old 5' v1 chemistry, `cellranger count` only, TCR/BCR libraries are in separate dataset pages (`vdj_v1_hs_pbmc2_t`, `_b`, UNVERIFIED). Not recommended.
- 10k Human Diseased PBMCs (ALL) Freshly Processed: https://www.10xgenomics.com/datasets/10k_5p_Human_diseased_PBMC_ALL_Fresh  (Cell Ranger 8.0.1, VDJ; page says only a subset of output files is provided; no FASTQs expected). UNVERIFIED, not suitable.

## Non-10x leads (GEO / ArrayExpress), not verified for raw FASTQ availability

- Stephenson et al. 2021, Nat Med, "Single-cell multi-omics analysis of the immune response in COVID-19": ArrayExpress **E-MTAB-10026**
  (https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-10026 resolves, title "Deciphering the molecular immune response to COVID-19 using single cell multi-omics").
  10x 5' GEX + CITE-seq + BCR + TCR from PBMCs of healthy and COVID-19 donors (from my memory of the paper: ~780k cells, ~190 Abs; those two numbers UNVERIFIED). UNVERIFIED also: whether raw FASTQs are open (access may be via EGA/controlled), which Cell Ranger version, and how libraries map to samples. Too large and too complex for v0.2 anyway.
- A search of GEO/SRA for "TotalSeq-C + TCR + BCR + PBMC + FASTQ" did not turn up an accession I could confirm in this session. If a non-10x dataset is wanted, searching GEO for `"cell surface protein" "VDJ" PBMC 10x 5'` with ENA fastq links is the next step. Not pursued further because the 10x sets above already satisfy every requirement with open FASTQs and published multi configs.

---

## Literature: B/T "dual expressers" versus doublets

1. Ahmed R, Omidian Z, Giwa A, Cornwell B, et al. *A public BCR present in a unique dual-receptor-expressing lymphocyte from type 1 diabetes patients encodes a potent T cell autoantigen.* **Cell 2019;177:1583-1599.e16**, doi:10.1016/j.cell.2019.05.007. Original claim: rare lymphocytes co-expressing TCR and BCR (dual expressers, DEs) enriched in T1D; scRNA-seq data GEO GSE129112 (Smart-seq2, 34 DEs).
2. Burel JG, Pomaznoy M, Lindestam Arlehamn CS, Seumois G, et al. *The challenge of distinguishing cell-cell complexes from singlet cells in non-imaging flow cytometry and single-cell sorting.* **Cytometry A 2020;97:1127-1135**, doi:10.1002/cyto.a.24027. Shows T cell-monocyte (and similar) complexes survive scatter gating and mimic dual-lineage cells.
3. Japp AS, Meng W, Rosenfeld AM, Perry DJ, et al. *TCR+/BCR+ dual-expressing cells and their associated public BCR clonotype are not enriched in type 1 diabetes.* **Cell 2021;184:827-839.e14**, doi:10.1016/j.cell.2020.11.035 (published online Dec 2020, issue Feb 2021). Matters Arising: DEs are not enriched in T1D, the public clonotype is not found in sorted cells. Ahmed et al. reply: Cell 2021;184:840-843, doi:10.1016/j.cell.2020.11.036.
4. Burel JG, Chawla S, Greenbaum JA, Peters B. *Distinguishing cell-cell complexes from dual lineage cells using single-cell transcriptomics is not trivial.* **Cytometry A 2022;101:547-551**, doi:10.1002/cyto.a.24656. Reanalysis of GSE129112: DE separation in PCA driven by high-gene-count outliers with doublet-like signatures. The defence is Jie C, Ahmed R, Hamad ARA, *Expression of unique gene signature distinguishes TCRab+/BCR+ dual expressers from CD3+CD14+ doublets*, Cytometry A 2022, doi:10.1002/cyto.a.24542.
5. Sun B, Bugarin-Estrada E, Overend LE, Walker CE, et al. *Double-jeopardy: scRNA-seq doublet/multiplet detection using multi-omic profiling.* **Cell Rep Methods 2021;1:100008**, doi:10.1016/j.crmeth.2021.100008. MLtiplet: uses co-capture of IG and TR (VDJ) or CITE-seq hybrids (e.g. CD19+CD3+) to flag doublets. The most directly relevant method for the planned check.
6. Recent (2024-2026), thinner evidence:
   - Bashford-Rogers, Shohan, Tucci et al., *CELLeBT: exploring physical B cell-T cell clonal interactions through single-cell multi-omics* (Research Square preprint, 18 Jun 2024, doi:10.21203/rs.3.rs-4437721/v1). Treats droplets with both BCR and TCR as candidate physical B-T interactions; preprint only, I did not find a journal version.
   - Yao ... Xu, *scRNA+TCR+BCR-seq revealed the proportions and gene expression patterns of dual receptor T and B lymphocytes in NPC and NLH*, **Biochem Biophys Res Commun 2024;709:149820**, doi:10.1016/j.bbrc.2024.149820. Reports dual-receptor lymphocytes as a biological population; relevance to doublet handling is limited.
   - I could not find a 2025-2026 paper that specifically re-tests B/T hybrids with a modern doublet detector (scDblFinder, Scrublet, DoubletFinder) against protein ground truth. Treat that as a gap this pipeline could fill. Note the task brief said "Burel et al. 2020 ... in Cell"; the Burel rebuttals I could verify are in Cytometry A (2020 and 2022), and Japp et al. is in Cell (e-pub 2020, vol 184 in 2021).

---

## Recommendation

**Best for v0.2: A, "Human PBMC from a Healthy Donor, 10k cells (v2)"** (10x, Cell Ranger 5.0.0, CC BY 4.0).

- Only candidate that is both small enough (**15.9 GB**, three tarballs) and complete: GEX 10,548 cells + Antibody Capture (CD3, CD19, CD4, CD8a, CD14, CD16, CD56, ... ) + VDJ-B (976 cells) + VDJ-T (5,308 cells).
- Best B-cell yield (976) of the practical options, so the B/T hybrid cross-tab has real counts: 95 barcodes already carry both IG and TR contigs in 10x's own filtered contigs.
- 10x publishes a ready `cellranger multi` config and feature reference (both OK 200). Same donor-type and 5' v2 chemistry/read layout as v0.1, with sequencing 26x90, so only minor config changes are needed for CR 10.0 (`reference` keys, `create-bam`, optionally the 2024-A / VDJ 7.1.0 references).
- Protein ground truth: CD3 for T, CD19 for B (no CD20; CD19 alone is adequate). Dual CD3+CD19 positive non-doublet cells (by RNA) are the hybrid test cases.

Config sketch (adapt paths; keys per CR 10 `multi`):

```
[gene-expression]
reference,/ref/refdata-gex-GRCh38-2024-A
create-bam,false
[vdj]
reference,/ref/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0
[feature]
reference,sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_feature_ref.csv
[libraries]
fastq_id,fastqs,lanes,feature_types
sc5p_v2_hs_PBMC_10k_5gex,<gex dir>,1|2,gene expression
sc5p_v2_hs_PBMC_10k_5fb,<5fb dir>,1|2,antibody capture
sc5p_v2_hs_PBMC_10k_b,<b dir>,1|2,vdj-b
sc5p_v2_hs_PBMC_10k_t,<t dir>,1|2,vdj-t
```

**Backup: B, "5k Human PBMCs stained with TotalSeq-C Human TBNK Cocktail, NextGEM 5'"** (10x, Cell Ranger 8.0.0, CC BY 4.0).

- 4,661 cells, GEX + Antibody (CD3, CD19, CD4, CD8, CD14, CD16, CD56, CD11c, CD45, IgG1) + VDJ-B (329) + VDJ-T (2,458), 10x config already in CR 8+ format with 2024-A and VDJ 7.1.0 (closest to CR 10).
- Costs: **33.1 GB** FASTQ download and fewer B cells. Choose this if a modern reference/Cell Ranger pairing matters more than download size, or if A's gex tar turns out to have an unexpected layout.

Not recommended for v0.2: C (74 GB), D (no protein), E (about 244 GB, but revisit later for hashtag-derived doublet truth and CD20), F (no TCR).
