# FUSION image for pheast's TWAS rules (twas_compute_weights_batch,
# twas_assemble_summary): the bundled FUSION R scripts and gcta_nr_robust from
# pheast/scripts/fusion_twas/ under /opt/fusion (utils/ layout kept, as
# FUSION.compute_weights.R sources utils/plink_utils.R relative to itself), and
# the conda env of pheast/envs/fusion.yaml (R 4.5.3, glmnet, optparse, Matrix,
# plink 1.90, gemma, Python 3.12), versions as in the baseline env
# (tests/hpc/envs/baseline-pheast.lock.txt).
# Built from the repository root (context .) by .github/workflows/fusion-image.yml
# and pushed to ghcr.io/twillis209/pantry-fusion.
# mambaorg/micromamba:2.10.0-ubuntu24.04
FROM mambaorg/micromamba@sha256:c3897a06f87c1c13f038188b3665a8c4162da62eac5109ae6a939d1e290dcefa

LABEL org.opencontainers.image.source="https://github.com/twillis209/Pantry" \
      org.opencontainers.image.description="FUSION TWAS scripts + gcta_nr_robust + R/plink/gemma for Pantry pheast" \
      org.opencontainers.image.licenses="MIT"

COPY --chown=$MAMBA_USER:$MAMBA_USER pheast/envs/fusion.yaml /tmp/fusion.yaml

RUN micromamba install -y -n base --override-channels -c conda-forge -c bioconda \
        -f /tmp/fusion.yaml \
    && micromamba clean -a -y

COPY --chown=$MAMBA_USER:$MAMBA_USER pheast/scripts/fusion_twas/ /opt/fusion/

# gcta_nr_robust is upstream Pantry 8da7787's binary (statically linked).
RUN echo "7b51e51d53d520589ce3c922cb13f4562cb67f842ffe6172aec7dfa1b91c93e9  /opt/fusion/gcta_nr_robust" \
        | sha256sum -c - \
    && chmod +x /opt/fusion/gcta_nr_robust

# apptainer exec skips the micromamba entrypoint, so put the env on PATH
# directly (python, Rscript, plink, gemma).
ENV PATH=/opt/conda/bin:$PATH

# Build-time gate: exact baseline versions, gcta links and runs, R packages
# load, plink and gemma run.
COPY containers/fusion/check_image.py /opt/check_image.py
RUN python /opt/check_image.py
