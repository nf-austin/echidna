# Single image for every process in this pipeline: the echidna stack is a
# superset of what SEG_TO_GENE_CN, CALL_TUMOR_CELLS and CONCAT_H5ADS need.
#
# Keep the pins in step with the modules' environment.yml files, which are the
# conda fallback for the same processes.
FROM python:3.10-slim

# wget + gzip for DOWNLOAD_GENE_BED, procps so Nextflow can collect task metrics.
RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        gzip \
        procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    "scanpy>=1.10" \
    "pyro-ppl>=1.9.1" \
    hmmlearn \
    leidenalg \
    python-igraph \
    seaborn \
    scikit-learn \
    pandas \
    anndata \
    PhenoGraph \
    "sc-echidna>=1.0.3"

# Fail the build rather than a run if an import is broken.
RUN python3 -c 'import scanpy, anndata, pandas, pyro, hmmlearn, sklearn; print("All packages OK")'
