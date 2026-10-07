#!/usr/bin/env bash
# One-shot runner for the GLOOM 0.2.0 revision on a server (no conda needed; Python 3.12 + requirements already installed).
# Usage (from anywhere):   PYTHON=/path/to/python bash scripts/run_remote.sh        # default: python
#   BENCH=1  also run the runtime benchmark at the end        SKIP_TESTS=1  skip pytest
# Safe to re-run: the GDC download is cached and the pipeline restarts from step 1.
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${PYTHON:-python}"
export PYTHONIOENCODING=utf-8

echo "== environment =="
$PY --version
$PY - <<'EOF'
import importlib, sys
need = ["numpy", "pandas", "scipy", "sklearn", "statsmodels", "xgboost", "networkx", "click", "requests", "gseapy", "matplotlib"]
bad = []
for m in need:
    try:
        mod = importlib.import_module(m); print(f"  {m:12s} {getattr(mod, '__version__', 'ok')}")
    except Exception as e:
        bad.append(m); print(f"  {m:12s} MISSING ({e})")
if bad:
    sys.exit("Missing packages: " + ", ".join(bad) + "  ->  pip install " + " ".join(bad))
EOF
python3 -c "import os; print('  cores:', os.cpu_count())" 2>/dev/null || true

if [ "${SKIP_TESTS:-0}" != "1" ]; then
  echo "== unit tests =="
  $PY -m pytest -q
fi

echo "== GDC TCGA-LUAD (tumor + adjacent normal, paired) =="
if [ ! -f data/raw/tcga_gdc/counts_matrix.csv ]; then
  $PY scripts/fetch_gdc_tcga_luad.py --paired-only --workers 6
else
  echo "  counts_matrix.csv already present, skipping download"
fi

echo "== config in use =="
grep -n '^DATA_SOURCE\|^DE_METHOD\|^SEED\|^USE_CROSSFIT\|^LABEL_INDEPENDENT_UNIVERSE\|^USE_BATCH_CORRECTION\|^CROSSFIT_REPEATS\|^BOOTSTRAP_N\|^ABLATION_REPEATS\|^NETWORK_BOOTSTRAP_N' src/gloom/pipeline/config.py

echo "== full pipeline =="
rm -rf outputs
$PY src/gloom/pipeline/run_pipeline.py 2>&1 | tee run_full.log

if [ "${BENCH:-0}" = "1" ]; then
  echo "== runtime benchmark =="
  $PY scripts/benchmark_runtime.py --genes 2000 5000 --samples 100 300 --repeats 3 || echo "benchmark failed (non-fatal)"
fi

echo "== return package =="
$PY -m pip freeze > pip_freeze.txt 2>/dev/null || true
R=return_package; rm -rf $R; mkdir -p $R/results $R/figures $R/data
cp run_full.log pip_freeze.txt $R/ 2>/dev/null || true
cp src/gloom/pipeline/config.py $R/config_used.py
cp outputs/logs/pipeline.log $R/ 2>/dev/null || true
cp -r outputs/results/*.csv outputs/results/*.txt $R/results/ 2>/dev/null || true
cp -r outputs/results/enrichment $R/results/ 2>/dev/null || true
cp -r outputs/results/reports $R/results/ 2>/dev/null || true
cp outputs/figures/*.png $R/figures/ 2>/dev/null || true
cp data/raw/tcga_gdc/sample_sheet.csv data/raw/tcga_gdc/gdc_files_manifest.csv $R/data/ 2>/dev/null || true
( command -v zip >/dev/null && zip -qr gloom_rerun_package.zip $R ) || tar -czf gloom_rerun_package.tar.gz $R
ls -la gloom_rerun_package.*
echo "DONE — send gloom_rerun_package.* back."
