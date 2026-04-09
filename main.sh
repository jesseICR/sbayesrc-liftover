#!/usr/bin/env bash
# =============================================================================
# SBayesRC snp.info hg19 -> hg38 Liftover
# =============================================================================
# Single entry point. Sets up a Python virtual environment in tools/venv/,
# installs dependencies, and runs the pipeline.
#
# Usage:
#   bash main.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/tools/venv"
PYTHON="$VENV_DIR/bin/python3"

# ---- Python virtual environment ---------------------------------------------
if [ -x "$PYTHON" ]; then
    echo "[skip] Python venv"
else
    echo "[setup] Creating Python venv in tools/venv/ ..."
    mkdir -p "$SCRIPT_DIR/tools"
    python3 -m venv "$VENV_DIR"
    "$VENV_DIR/bin/pip" install --upgrade pip --quiet
    "$VENV_DIR/bin/pip" install --no-cache-dir -r "$SCRIPT_DIR/requirements.txt" --quiet
    echo "[done] Python venv ready"
fi

# ---- Logging -- all subsequent output goes to both terminal and log file ----
mkdir -p "$SCRIPT_DIR/logs"
LOG_FILE="$SCRIPT_DIR/logs/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

# ---- Run pipeline ------------------------------------------------------------
"$PYTHON" -u "$SCRIPT_DIR/liftover.py" "$@"
