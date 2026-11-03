#!/usr/bin/env python
"""Check that every target of one dataset exists and is non-empty.

Usage: check_targets.py <phenotyping|pheast> <config.yml> <dataset> <output_dir>

Targets come from the stage's own lib/config.py (normalise + dataset_targets),
so the list is the one 'rule all' requests. Per target: exists and size > 0;
*.bed.gz, *.txt.gz (QTL tables), *.tsv: header plus >= 1 data row;
*.tar.bz2: >= 1 regular file; other text files (phenotype groups): >= 1
line. Prints one line per target and exits 1 if any fails.
"""
import gzip
import importlib.util
import sys
import tarfile
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]


def load_lib(stage: str):
    spec = importlib.util.spec_from_file_location(
        f"{stage}_config", REPO / stage / "lib" / "config.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def n_lines(path: Path) -> int:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as f:
        return sum(1 for _ in f)


def check(path: Path) -> tuple:
    """Return (ok, description) for one target."""
    if not path.exists():
        return False, "missing"
    size = path.stat().st_size
    if size == 0:
        return False, "empty file"
    name = path.name
    if name.endswith(".tbi"):
        return True, f"{size} bytes"
    if name.endswith(".tar.bz2"):
        with tarfile.open(path, "r:bz2") as tar:
            files = [m for m in tar.getmembers() if m.isfile()]
        return len(files) > 0, f"{len(files)} files in archive"
    lines = n_lines(path)
    if name.endswith((".bed.gz", ".txt.gz", ".tsv")):
        return lines > 1, f"{lines - 1} data rows"
    return lines > 0, f"{lines} lines"


def main() -> None:
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    stage, config_path, dataset, output_dir = sys.argv[1:]
    lib = load_lib(stage)
    config = lib.normalise(yaml.safe_load(Path(config_path).read_text()))
    targets = lib.dataset_targets(config["datasets"][dataset])
    failed = 0
    for target in targets:
        ok, what = check(Path(output_dir) / target)
        failed += not ok
        print(f"{'OK  ' if ok else 'FAIL'} {target}: {what}")
    print(f"{stage} {dataset}: {len(targets) - failed}/{len(targets)} targets OK")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
