"""GRCh38 ref/alt and Open Targets IDs for a genotype bim (R36).

Writes variant_id, chrom, pos, ref, alt, ot_variant_id, flip and status
(susie_tables.orient_alleles), one chromosome at a time, so memory is about
one chromosome's sequence. The FASTA is streamed: no index is read or written.
  variant_alleles.py BIM FASTA OUT COUNTED
COUNTED is a1 (bim column 5) or a2 (column 6), the allele tensorQTL counts.
"""

import sys

import pandas as pd

import susie_tables  # beside this script


def main(bim_path, fasta, out, counted):
    bim = pd.read_csv(
        bim_path, sep=r"\s+", header=None, usecols=[0, 1, 3, 4, 5],
        names=["chrom", "snp", "pos", "a1", "a2"], dtype={0: str, 1: str, 4: str, 5: str},
    )
    bim["key"] = bim["chrom"].str.removeprefix("chr")
    by_chrom = dict(tuple(bim.groupby("key", sort=False)))
    parts = []
    for name, seq in susie_tables.read_fasta(fasta, keep=set(by_chrom)):
        def ref_base(chrom, pos, length, seq=seq):
            return seq[pos - 1 : pos - 1 + length]
        parts.append(susie_tables.orient_alleles(by_chrom.pop(name), ref_base, counted))
    if by_chrom:
        raise SystemExit(f"chromosomes not in {fasta}: {', '.join(sorted(by_chrom))}")
    table = pd.concat(parts)
    table = table.loc[bim["snp"]]  # bim order
    table.to_csv(out, sep="\t")
    counts = table["status"].value_counts()
    print("  * alleles: " + ", ".join(f"{k} {counts.get(k, 0)}" for k in ["ok", "flipped", "mismatch"]))


if __name__ == "__main__":
    main(*sys.argv[1:])
