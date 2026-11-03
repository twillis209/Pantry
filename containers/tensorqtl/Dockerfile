# tensorQTL 1.0.10 / torch 2.3.1 (CUDA 12.1 wheel) image for pheast's tensorQTL
# rules, matching the baseline env (tests/hpc/envs/baseline-pheast.lock.txt).
# Built and pushed by .github/workflows/tensorqtl-image.yml to
# ghcr.io/twillis209/pantry-tensorqtl. The torch wheel bundles the CUDA
# runtime; the host driver comes in via apptainer --nv.
# mambaorg/micromamba:2.10.0-ubuntu24.04
FROM mambaorg/micromamba@sha256:c3897a06f87c1c13f038188b3665a8c4162da62eac5109ae6a939d1e290dcefa

LABEL org.opencontainers.image.source="https://github.com/twillis209/Pantry" \
      org.opencontainers.image.description="tensorQTL 1.0.10 + torch 2.3.1 (cu121) for Pantry pheast" \
      org.opencontainers.image.licenses="MIT"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml constraints.txt /tmp/

RUN micromamba install -y -n base --override-channels -c conda-forge -c bioconda \
        -f /tmp/environment.yml \
    && /opt/conda/bin/python -m pip install --no-cache-dir -c /tmp/constraints.txt \
        tensorqtl==1.0.10 torch==2.3.1 \
    && micromamba clean -a -y

# apptainer exec skips the micromamba entrypoint, so put the env on PATH
# directly (python3, R for rpy2) and name R_HOME for rpy2.
ENV PATH=/opt/conda/bin:$PATH \
    R_HOME=/opt/conda/lib/R

# Build-time gate: exact baseline versions, torch<->numpy 2 interop, qvalue.
COPY check_image.py /opt/check_image.py
RUN python3 /opt/check_image.py
