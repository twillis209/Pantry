# MAJIQ v3 image for Pantry's intron_retention rules (majiq_container).
#
# LICENCE: MAJIQ is under the University of Pennsylvania academic licence,
# https://majiq.biociphers.org/Academic_License.txt. You must accept that
# licence before building this image. Use is for non-commercial or educational
# purposes only. Do NOT redistribute the built image (do not push it to a
# public or shared registry, or share it): each user builds their own. This
# recipe contains no MAJIQ code; MAJIQ is fetched at build time from
# https://bitbucket.org/biociphers/majiq_academic (tag v3.0.22, commit d841baad).
# No Pantry workflow builds or publishes this image.
#
# For machines with Docker but no apptainer; majiq.def (apptainer) runs the
# same steps. Build and convert to a .sif for the HPC:
#   docker build --platform linux/amd64 -t pantry-majiq:3.0.22 containers/majiq
#   docker save pantry-majiq:3.0.22 -o majiq.tar
#   apptainer build majiq.sif docker-archive://majiq.tar
# mambaorg/micromamba:2.10.0-ubuntu24.04, as the other Pantry images.
FROM mambaorg/micromamba@sha256:c3897a06f87c1c13f038188b3665a8c4162da62eac5109ae6a939d1e290dcefa

LABEL org.opencontainers.image.title="pantry-majiq (user-built, not for redistribution)" \
      org.opencontainers.image.source="https://github.com/twillis209/Pantry" \
      org.opencontainers.image.version="3.0.22" \
      org.opencontainers.image.licenses="LicenseRef-Penn-MAJIQ-Academic" \
      majiq.source="https://bitbucket.org/biociphers/majiq_academic" \
      majiq.commit="d841baad8a3498251e0113fe5f7c41783c4eb55f" \
      majiq.licence="https://majiq.biociphers.org/Academic_License.txt" \
      majiq.copyright="Copyright © 2016 University of Pennsylvania; All Rights Reserved" \
      majiq.creators="Yoseph Barash, Jorge V. Garcia, Alejandro Barrera, Matthew Gazzara" \
      majiq.notices="/opt/majiq"

COPY --chown=$MAMBA_USER:$MAMBA_USER majiq3.yml majiq3-install.sh /tmp/majiq-build/

# The env of the conda route (majiq3.yml), as the image's base env.
RUN micromamba install -y -n base --override-channels -c conda-forge -c bioconda \
        -f /tmp/majiq-build/majiq3.yml

USER root
RUN mkdir -p /opt/majiq && chown $MAMBA_USER:$MAMBA_USER /opt/majiq
USER $MAMBA_USER

# The same install steps as the conda route; notices to /opt/majiq.
RUN bash /tmp/majiq-build/majiq3-install.sh --prefix /opt/conda --notices /opt/majiq \
        /tmp/majiq-build/work \
    && rm -rf /tmp/majiq-build \
    && micromamba clean -a -y

# apptainer exec does not run the micromamba entrypoint: put the env on PATH.
ENV PATH=/opt/conda/bin:$PATH

RUN majiq --version | grep -q 3.0.22 && python -c "import rna_majiq, rna_voila"
