#!/usr/bin/env bash
# MAJIQ v3 is under Penn's academic licence
# (https://majiq.biociphers.org/Academic_License.txt). Accept it before running
# this script. Use is for non-commercial or educational purposes only. Do not
# redistribute MAJIQ, or any env or image built with it.
#
# Install MAJIQ v3 (academic, https://bitbucket.org/biociphers/majiq_academic,
# tag v3.0.22) by "pip install ." from the repository root, with HTSLIB_* set,
# as in its README. This is the single set of install steps for both routes:
#
#   bash majiq3-install.sh <workdir>
#       conda route: (re)create the named conda env 'majiq3' from majiq3.yml,
#       then install MAJIQ into it (Pantry's majiq_env default).
#   bash majiq3-install.sh --prefix <env> <workdir>
#       image builds (majiq.def, Dockerfile): install into the existing env at
#       <env>, already created from majiq3.yml.
#
# MAJIQ is not patched. The build stops unless the tag is at the expected
# commit. MAJIQ's copyright and licence notices go to
# <env>/share/majiq (--notices <dir> to change).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
tag=v3.0.22
commit=d841baad8a3498251e0113fe5f7c41783c4eb55f
version=3.0.22
url=https://bitbucket.org/biociphers/majiq_academic.git

prefix=
notices=
while [ $# -gt 1 ]; do
    case $1 in
        --prefix) prefix=$2; shift 2 ;;
        --notices) notices=$2; shift 2 ;;
        *) break ;;
    esac
done
work=${1:?usage: majiq3-install.sh [--prefix <env>] [--notices <dir>] <workdir>}
mkdir -p "$work"
work=$(cd "$work" && pwd)

if [ -z "$prefix" ]; then
    mamba env remove -y -n majiq3 >/dev/null 2>&1 || true
    mamba env create -y -f "$here/majiq3.yml" --override-channels -c conda-forge -c bioconda
    eval "$(conda shell.bash hook)"
    conda activate majiq3
else
    # Activate (not just PATH) so the conda compilers set CC/CXX.
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate "$prefix"
fi
notices=${notices:-$CONDA_PREFIX/share/majiq}
export HTSLIB_LIBRARY_DIR=$CONDA_PREFIX/lib HTSLIB_INCLUDE_DIR=$CONDA_PREFIX/include
# Let CMake find the env's libraries before the system ones.
export CMAKE_PREFIX_PATH=$CONDA_PREFIX

src=$work/majiq_academic
rm -rf "$src"
# Full clone: setuptools_scm derives the version from the tag.
git clone -q "$url" "$src"
git -C "$src" checkout -q "$tag"
head=$(git -C "$src" rev-parse HEAD)
echo "majiq_academic $tag commit $head"
if [ "$head" != "$commit" ]; then
    echo "error: tag $tag is at $head, expected $commit" >&2
    exit 1
fi
(cd "$src" && pip install --no-cache-dir .)

# Keep MAJIQ's notices (licence Section 7): every licence/copying file of the
# clone, at its path in the repository, plus Penn's copyright notice.
mkdir -p "$notices"
(cd "$src" && find . -path ./.git -prune -o -type f \
    \( -iname 'LICENSE*' -o -iname 'COPYING*' -o -iname 'NOTICE*' \) \
    -exec cp --parents {} "$notices" \;)
cat >"$notices/NOTICE" <<EOF
MAJIQ and VOILA $version ($url, tag $tag, commit $commit)
Copyright © 2016 University of Pennsylvania; All Rights Reserved
Creators: Yoseph Barash, Jorge V. Garcia, Alejandro Barrera, Matthew Gazzara
(School of Medicine, University of Pennsylvania).
Used under the Penn academic licence, https://majiq.biociphers.org/Academic_License.txt:
non-commercial or educational use only; do not redistribute.
The other files here are the licences of the third-party code bundled with MAJIQ.
EOF

# Build check: the installed version and the Python packages.
got=$(majiq --version 2>&1 | tail -n 1)
echo "majiq --version: $got"
case $got in
    *"$version"*) ;;
    *) echo "error: majiq --version is '$got', expected $version" >&2; exit 1 ;;
esac
python -c "import rna_majiq, rna_voila; print('import rna_majiq, rna_voila: ok')"
