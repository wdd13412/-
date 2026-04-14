#!/usr/bin/env bash
set -euo pipefail

echo "[env-setup] Updating apt index..."
sudo apt-get update -qq

echo "[env-setup] Installing build dependencies (gfortran/lapack/blas)..."
sudo apt-get install -y -qq gfortran liblapack-dev libblas-dev

echo "[env-setup] Toolchain ready:"
gfortran --version | head -n 1
