#!/usr/bin/env bash
# Create the named conda env 'majiq3' and pip-install MAJIQ v3 (academic,
# https://bitbucket.org/biociphers/majiq_academic) into it, following the
# repository README ("pip install ." from the repo root, HTSLIB_* set).
# Usage (inside an sjob): bash majiq3-install.sh <workdir>
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
work=${1:?workdir}
tag=v3.0.22
mkdir -p "$work"
mamba env remove -y -n majiq3 >/dev/null 2>&1 || true
mamba env create -y -f "$here/majiq3.yml" --override-channels -c conda-forge -c bioconda
eval "$(conda shell.bash hook)"
conda activate majiq3
export HTSLIB_LIBRARY_DIR=$CONDA_PREFIX/lib HTSLIB_INCLUDE_DIR=$CONDA_PREFIX/include
# Let CMake find the env's libraries before the system ones.
export CMAKE_PREFIX_PATH=$CONDA_PREFIX
rm -rf "$work/majiq_academic"
# Full clone: setuptools_scm derives the version from the tag.
git clone -q https://bitbucket.org/biociphers/majiq_academic.git "$work/majiq_academic"
git -C "$work/majiq_academic" checkout -q "$tag"
echo "majiq_academic $tag commit $(git -C "$work/majiq_academic" rev-parse HEAD)"
cd "$work/majiq_academic"
pip install --no-cache-dir .
cd "$work"
majiq --version
