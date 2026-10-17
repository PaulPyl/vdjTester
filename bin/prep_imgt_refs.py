#!/usr/bin/env python3
"""Derive every tool's germline input from one IMGT snapshot (FETCH_IMGT).

The snapshot has two parts, downloaded together:
  bulk.fa                  IMGT/GENE-DB reference sequences, F + ORF + in-frame P,
                           IMGT-gapped V/D/J regions (all species; human is used)
  mkvdjref-input.fasta     Cell Ranger's fetch-imgt output: leader + V (L-PART1 +
                           V-EXON), D, J and C, in `cellranger mkvdjref --seqs` format

Cell Ranger needs the leader to call a contig full-length, and IMGT has no leader
for some alleles; the bulk set leaves out pseudogenes. So all tools get the same V,
D and J alleles: those in both parts. On IMGT of 2026-10-10 that is 1,098 alleles of
511 genes (bulk: 1,329 of 523; the 12 genes lost are mostly orphons such as
IGHV1/OR15-1). Constant genes are not compared between tools, so each keeps its
usual source: IMGT C for Cell Ranger, MiXCR's built-in C genes, TRUST4's genomic C.

Outputs in --outdir:
  cellranger/seqs.fasta                       mkvdjref --seqs input
  germlines/human/vdj/imgt_human_<X>.fasta    IMGT-gapped, for imgt2igblast / Change-O
  mixcr/<locus><segment>.fasta                ungapped, for mixcr buildLibrary
  trust4/human_IMGT+C.fa                      gapped V/D/J + TRUST4's C genes (--ref)
  igblast/human_gl.aux                        IgBLAST J frames and CDR3 ends for these J
                                              alleles (the container's aux file lacks newer
                                              ones, e.g. IGLJ2A/2B, and IgBLAST then reports
                                              no junction for them)
  alleles.tsv                                 every allele, which part has it, used or not
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

VDJ_LABELS = {"V-REGION": "V", "D-REGION": "D", "J-REGION": "J"}
CR_VDJ_TYPES = {"L-REGION+V-REGION", "D-REGION", "J-REGION"}


def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def read_bulk(path):
    """Human IG/TR V/D/J records: allele -> (header fields, gapped sequence).
    Duplicate records of an allele (IGLV3-21*05 appears twice) keep the first."""
    alleles = {}
    for header, seq in read_fasta(path):
        f = header.split("|")
        if len(f) < 5 or f[2] != "Homo sapiens" or f[4] not in VDJ_LABELS or f[1][:2] not in ("IG", "TR"):
            continue
        allele = f[1]
        if allele in alleles:
            if alleles[allele][1].replace(".", "").upper() != seq.replace(".", "").upper():
                print(f"{allele}: duplicate records differ, kept the first", file=sys.stderr)
            continue
        alleles[allele] = (f, seq.upper())
    return alleles


def read_cellranger(path):
    """fetch-imgt records: list of (header, seq, allele, region type)."""
    recs = []
    for header, seq in read_fasta(path):
        f = header.split("|")
        recs.append((header, seq, f[1].split()[0], f[3]))
    return recs


CODONS = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
J_CHAIN = {"IGH": "JH", "IGK": "JK", "IGL": "JL", "TRA": "JA", "TRB": "JB", "TRD": "JD", "TRG": "JG"}


def igblast_aux_line(allele, fields, seq):
    """IgBLAST aux row: allele, coding frame start, chain type, CDR3 stop, extra bases
    (0-based). The frame is IMGT's codon start; CDR3 ends just before the W/F of the
    conserved W/F-G-x-G motif. Without the motif (e.g. some P alleles) only the frame
    and chain type are given, as in IgBLAST's own file."""
    seq = seq.replace(".", "").upper()
    frame = int(fields[7]) - 1
    aa = "".join(CODONS.get(seq[i:i + 3], "X") for i in range(frame, len(seq) - 2, 3))
    chain = J_CHAIN[allele[:3]]
    for i in range(len(aa) - 3):
        if aa[i] in "WF" and aa[i + 1] == "G" and aa[i + 3] == "G":
            return f"{allele}\t{frame}\t{chain}\t{frame + 3 * i - 1}\t{(len(seq) - frame) % 3}"
    return f"{allele}\t{frame}\t{chain}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bulk", required=True)
    ap.add_argument("--cellranger-input", required=True)
    ap.add_argument("--trust4-default", required=True, help="TRUST4's human_IMGT+C.fa, for its C genes")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    out = Path(args.outdir)

    bulk = read_bulk(args.bulk)
    cr = read_cellranger(args.cellranger_input)
    cr_vdj = {allele for _h, _s, allele, rtype in cr if rtype in CR_VDJ_TYPES}
    common = set(bulk) & cr_vdj

    # alleles.tsv
    with open(out / "alleles.tsv", "w") as fh:
        fh.write("allele\tgene\tfunctionality\tin_imgt_bulk\tin_fetch_imgt\tused\n")
        for allele in sorted(set(bulk) | cr_vdj):
            func = bulk[allele][0][3] if allele in bulk else ""
            fh.write(f"{allele}\t{allele.split('*')[0]}\t{func}\t{allele in bulk}\t{allele in cr_vdj}\t{allele in common}\n")

    # Cell Ranger: shared V/D/J alleles, all C records.
    (out / "cellranger").mkdir(parents=True, exist_ok=True)
    with open(out / "cellranger" / "seqs.fasta", "w") as fh:
        for header, seq, allele, rtype in cr:
            if rtype not in CR_VDJ_TYPES or allele in common:
                fh.write(f">{header}\n{seq}\n")

    # Per locus and segment, from the bulk (gapped) sequences.
    groups = defaultdict(list)
    for allele in sorted(common):
        f, seq = bulk[allele]
        groups[(allele[:3], VDJ_LABELS[f[4]])].append((allele, f, seq))

    vdj = out / "germlines" / "human" / "vdj"
    vdj.mkdir(parents=True, exist_ok=True)
    (out / "mixcr").mkdir(exist_ok=True)
    for (locus, seg), recs in groups.items():
        # Change-O / imgt2igblast: IMGT headers, species with an underscore, gaps kept.
        with open(vdj / f"imgt_human_{locus}{seg}.fasta", "w") as fh:
            for _allele, f, seq in recs:
                f = list(f)
                f[2] = f[2].replace(" ", "_")
                fh.write(f">{'|'.join(f)}\n{seq}\n")
        # MiXCR: allele names, no gaps.
        with open(out / "mixcr" / f"{locus}{seg}.fasta", "w") as fh:
            for allele, _f, seq in recs:
                fh.write(f">{allele}\n{seq.replace('.', '')}\n")

    # TRUST4: gapped V/D/J (it reads IMGT positions from the gaps), then its own C genes.
    (out / "trust4").mkdir(exist_ok=True)
    with open(out / "trust4" / "human_IMGT+C.fa", "w") as fh:
        for allele in sorted(common):
            fh.write(f">{allele}\n{bulk[allele][1]}\n")
        n_c = 0
        for header, seq in read_fasta(args.trust4_default):
            name = header.split()[0].split("*")[0]
            if len(name) >= 4 and name[:2] in ("IG", "TR") and name[3] == "C":
                fh.write(f">{header}\n{seq}\n")
                n_c += 1

    (out / "igblast").mkdir(exist_ok=True)
    with open(out / "igblast" / "human_gl.aux", "w") as fh:
        fh.write("#gene/allele name, first coding frame start position, chain type, CDR3 stop, "
                 "extra bps beyond J coding end.\n#All positions are 0-based\n\n")
        for allele in sorted(a for a in common if VDJ_LABELS[bulk[a][0][4]] == "J"):
            fh.write(igblast_aux_line(allele, *bulk[allele]) + "\n")

    genes = {a.split("*")[0] for a in common}
    print(f"IMGT bulk: {len(bulk)} V/D/J alleles; fetch-imgt: {len(cr_vdj)}; shared: {len(common)} alleles of "
          f"{len(genes)} genes; TRUST4 C records: {n_c}")


if __name__ == "__main__":
    main()
