"""Synthetic configs for dry-run scaling measurements (Task 18).

Usage: make_scaling_config.py <rundir> <n_datasets> <n_samples>

Writes, under <rundir> (which needs the input/{phenotyping,pheast} links of
tests/host/config.yaml):
- fastq/S0001_{1,2}.fastq.gz ... : empty files (a dry run only needs them to exist);
- config.yaml: tests/host/config.yaml with phenotyping datasets d01..dNN, each
  a copy of 'full' with an inline sample list and FASTQ map, and pheast
  datasets d01..dNN chained to them (phenotype_dir results/pantry/phenotypes/dNN);
- standalone/: an old-format phenotyping config.yml (tests/configs/
  phenotyping_baseline.yml) with samples.txt and fastq_map.txt for the same
  samples, plus an input link, for one-dataset standalone runs (also upstream's).
"""

import copy
import os
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]


def main(rundir: str, n_datasets: int, n_samples: int):
    run = Path(rundir).resolve()
    fq = run / "fastq"
    fq.mkdir(parents=True, exist_ok=True)
    samples = [f"S{i:04d}" for i in range(1, n_samples + 1)]
    for s in samples:
        for r in (1, 2):
            (fq / f"{s}_{r}.fastq.gz").touch()

    host = yaml.safe_load(open(REPO / "tests" / "host" / "config.yaml"))
    pheno, pheast = host["pantry"]["phenotyping"], host["pantry"]["pheast"]
    template = pheno["datasets"]["full"]
    fastq_map = {s: [[f"{s}_1.fastq.gz"], [f"{s}_2.fastq.gz"]] for s in samples}
    names = [f"d{i:02d}" for i in range(1, n_datasets + 1)]
    pheno["datasets"] = {}
    pheast["datasets"] = {}
    for name in names:
        ds = copy.deepcopy(template)
        ds.update(samples=list(samples), fastq_map=fastq_map, fastq_dir=str(fq))
        pheno["datasets"][name] = ds
        pheast["datasets"][name] = {
            "phenotype_dir": f"results/pantry/phenotypes/{name}",
            "samples": list(samples),
        }
    with open(run / "config.yaml", "w") as f:
        yaml.safe_dump(host, f, sort_keys=False)

    sa = run / "standalone"
    sa.mkdir(exist_ok=True)
    legacy = yaml.safe_load(open(REPO / "tests" / "configs" / "phenotyping_baseline.yml"))
    legacy.update(samples_file="samples.txt", fastq_map="fastq_map.txt", fastq_dir=str(fq))
    with open(sa / "config.yml", "w") as f:
        yaml.safe_dump(legacy, f, sort_keys=False)
    (sa / "samples.txt").write_text("".join(f"{s}\n" for s in samples))
    (sa / "fastq_map.txt").write_text(
        "".join(f"{s}_1.fastq.gz\t{s}_2.fastq.gz\t{s}\n" for s in samples)
    )
    link = sa / "input"
    if not link.exists():
        os.symlink(os.readlink(run / "input" / "phenotyping"), link)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
