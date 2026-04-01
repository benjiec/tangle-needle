FROM mambaorg/micromamba:latest

USER root
RUN apt-get update && apt-get install -y \
    git

USER $MAMBA_USER
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    hmmer=3.4 \
    python=3.13 \
    pip=26.0.1 \
    && micromamba clean --all --yes

WORKDIR /needle

COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml ./ 
RUN micromamba run -n base pip install --no-cache-dir . || true

COPY --chown=$MAMBA_USER:$MAMBA_USER . .
RUN micromamba run -n base pip install --no-cache-dir .

ENV PATH="/needle:${PATH}"
ENV PYTHONPATH="/needle"

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# We don't use ENTRYPOINT or CMD because Nextflow 
# will override them to run its own wrapper script.
