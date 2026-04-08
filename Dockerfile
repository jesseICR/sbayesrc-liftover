# =============================================================================
# SBayesRC snp.info hg19 to hg38 Liftover — Docker Image
# =============================================================================
# Builds a minimal image with Python dependencies and curl pre-installed.
# Reference data (~5 GB) is downloaded at runtime and cached in the working
# directory.  Mount the repo directory so downloads persist across runs:
#
# Build:
#   docker build -t sbayesrc-liftover .
#
# Run:
#   docker run --rm -v $(pwd):/data sbayesrc-liftover
#
# Pull from GHCR instead of building locally:
#   docker pull ghcr.io/jesseicr/sbayesrc-liftover:latest
#   docker run --rm -v $(pwd):/data ghcr.io/jesseicr/sbayesrc-liftover:latest
#
# =============================================================================

FROM python:3.11-slim-bookworm

LABEL org.opencontainers.image.source="https://github.com/jesseICR/sbayesrc-liftover"
LABEL org.opencontainers.image.description="SBayesRC snp.info hg19 to hg38 liftover with dbSNP position verification"
LABEL org.opencontainers.image.licenses="MIT"

# ---- System dependencies ----------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# ---- Python dependencies (cached layer) -------------------------------------
COPY requirements.txt /app/
RUN pip install --no-cache-dir -r /app/requirements.txt

# ---- Pipeline code -----------------------------------------------------------
COPY liftover.py /app/

WORKDIR /data
CMD ["python3", "-u", "/app/liftover.py"]
