#!/usr/bin/env python
"""Write the derived inputs for the multi-dataset test configs.

Usage:
  make_multi_inputs.py phenotyping <input_dir> <out_dir> [--n 8]
  make_multi_inputs.py pheast <input_dir> <out_dir> [--n 300]

phenotyping: <out_dir>/samples.txt holds the first N samples of
<input_dir>/samples.txt (file order) and <out_dir>/fastq_map.txt their rows of
<input_dir>/fastq_map.txt with R1 only ('R1<TAB>sample', paths still relative
to the dataset's fastq_dir), i.e. the paired-end test data read as single-end.

pheast: <out_dir>/samples.txt holds the first N samples of
<input_dir>/samples.txt; <out_dir>/phenotypes/ has every
<input_dir>/phenotypes/*.bed.gz cut to its first four columns plus those
samples (in BED column order), bgzipped and tabix-indexed, and a copy of every
*.phenotype_groups.txt. pheast's samples list does not subset the analysis, so
the subset dataset needs BEDs with only its samples.

Deterministic: no sorting, no randomness. Needs bgzip and tabix on PATH.
"""
import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path


def first_samples(path: Path, n: int) -> list:
    samples = [line for line in path.read_text().splitlines() if line.strip()]
    if len(samples) < n:
        sys.exit(f"{path} has {len(samples)} samples, fewer than {n}")
    return samples[:n]


def write_lines(path: Path, lines: list) -> None:
    path.write_text("".join(f"{line}\n" for line in lines))


def phenotyping(input_dir: Path, out_dir: Path, n: int) -> None:
    samples = first_samples(input_dir / "samples.txt", n)
    keep = set(samples)
    rows = []
    for line in (input_dir / "fastq_map.txt").read_text().splitlines():
        fields = line.split("\t")
        if fields[-1] in keep:
            rows.append(f"{fields[0]}\t{fields[-1]}")
    missing = keep - {row.split("\t")[-1] for row in rows}
    if missing:
        sys.exit(f"samples without FASTQs: {', '.join(sorted(missing))}")
    out_dir.mkdir(parents=True, exist_ok=True)
    write_lines(out_dir / "samples.txt", samples)
    write_lines(out_dir / "fastq_map.txt", rows)
    print(f"phenotyping: {len(samples)} samples, {len(rows)} single-end FASTQ rows")


def subset_bed(src: Path, dest: Path, samples: list) -> int:
    with gzip.open(src, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        missing = set(samples) - set(header[4:])
        if missing:
            sys.exit(f"{src}: samples not in BED: {', '.join(sorted(missing))}")
        keep = set(samples)
        cols = list(range(4)) + [i for i in range(4, len(header)) if header[i] in keep]
        plain = dest.with_suffix("")  # x.bed.gz -> x.bed
        n_rows = 0
        with open(plain, "w") as out:
            out.write("\t".join(header[i] for i in cols) + "\n")
            for line in f:
                fields = line.rstrip("\n").split("\t")
                out.write("\t".join(fields[i] for i in cols) + "\n")
                n_rows += 1
    subprocess.run(["bgzip", "-f", str(plain)], check=True)
    subprocess.run(["tabix", "-f", "-p", "bed", str(dest)], check=True)
    return n_rows


def pheast(input_dir: Path, out_dir: Path, n: int) -> None:
    samples = first_samples(input_dir / "samples.txt", n)
    pheno_in = input_dir / "phenotypes"
    pheno_out = out_dir / "phenotypes"
    pheno_out.mkdir(parents=True, exist_ok=True)
    write_lines(out_dir / "samples.txt", samples)
    beds = sorted(pheno_in.glob("*.bed.gz"))
    if not beds:
        sys.exit(f"no *.bed.gz in {pheno_in}")
    for bed in beds:
        n_rows = subset_bed(bed, pheno_out / bed.name, samples)
        print(f"pheast: {bed.name}: {n_rows} phenotypes x {len(samples)} samples")
    for groups in sorted(pheno_in.glob("*.phenotype_groups.txt")):
        shutil.copyfile(groups, pheno_out / groups.name)
        print(f"pheast: copied {groups.name}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("stage", choices=["phenotyping", "pheast"])
    ap.add_argument("input_dir", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--n", type=int, help="number of samples (default 8 / 300)")
    args = ap.parse_args()
    for tool in ("bgzip", "tabix") if args.stage == "pheast" else ():
        if shutil.which(tool) is None:
            sys.exit(f"{tool} not found on PATH")
    if args.stage == "phenotyping":
        phenotyping(args.input_dir, args.out_dir, args.n or 8)
    else:
        pheast(args.input_dir, args.out_dir, args.n or 300)


if __name__ == "__main__":
    main()
