"""Check the FUSION image: baseline versions, gcta, R packages, plink, gemma.

Run at build time (Dockerfile) and on the HPC.
"""

import glob
import json
import subprocess

GCTA = "/opt/fusion/gcta_nr_robust"


def run(*cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


# Conda package versions, as in tests/hpc/envs/baseline-pheast.lock.txt, from
# the env's package records.
versions = {}
for path in glob.glob("/opt/conda/conda-meta/*.json"):
    with open(path) as f:
        record = json.load(f)
    versions[record["name"]] = record["version"]
want = {
    "python": "3.12.16",
    "r-base": "4.5.3",
    "r-glmnet": "5.0",
    "r-optparse": "1.8.2",
    "r-matrix": "1.7_6",
    "plink": "1.90",
    "gemma": "0.98.5",
}
got = {name: versions.get(name) for name in want}
print(got)
assert got == want, f"versions differ from the baseline: {got} != {want}"

# gcta_nr_robust is statically linked: ldd reports "not a dynamic executable";
# a dynamic build must resolve every library.
ldd = run("ldd", GCTA)
print("ldd:", (ldd.stdout + ldd.stderr).strip())
assert "not found" not in ldd.stdout + ldd.stderr, "gcta has unresolved libraries"
assert ldd.returncode == 0 or "not a dynamic executable" in ldd.stdout + ldd.stderr
# FUSION's own check: system(PATH_gcta) must return 0.
assert run(GCTA).returncode == 0, "gcta_nr_robust does not run"

r = run("Rscript", "-e", "library(glmnet); library(optparse); library(Matrix)")
assert r.returncode == 0, f"R packages do not load: {r.stderr}"
print("R glmnet, optparse, Matrix: ok")

# FUSION's checks: plink --help and gemma -h must return 0.
assert run("plink", "--help").returncode == 0, "plink does not run"
assert run("gemma", "-h").returncode == 0, "gemma does not run"
print("gcta, plink, gemma: ok")
