#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

RUN_AFTER_BUILD=false
if [[ "${1:-}" == "--run" ]]; then
  RUN_AFTER_BUILD=true
fi

required_files=(
  "TypesModule.f90"
  "meshdeformationn_d.f90"
  "BuFlow_test_d.f90"
  "main_d.f90"
  "run_parameter_d_AAA.f90"
  "mesh/OFairfoilMesh/points"
  "mesh/OFairfoilMesh/faces"
  "mesh/OFairfoilMesh/owner"
  "mesh/OFairfoilMesh/neighbour"
  "mesh/OFairfoilMesh/boundary"
)

for f in "${required_files[@]}"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] Missing required file: $f" >&2
    exit 1
  fi
done

# Prefer gfortran-11 (user-verified working version), then fall back to FC/gfortran.
FC="${FC:-}"
if [[ -z "$FC" ]]; then
  if command -v gfortran-11 >/dev/null 2>&1; then
    FC="gfortran-11"
  elif command -v gfortran >/dev/null 2>&1; then
    FC="gfortran"
  fi
fi

if [[ -z "$FC" ]] || ! command -v "$FC" >/dev/null 2>&1; then
  echo "[ERROR] Fortran compiler not found (gfortran/gfortran-11)." >&2
  echo "Install dependencies first (Ubuntu/Debian):" >&2
  echo "  sudo apt update" >&2
  echo "  sudo apt install -y gfortran-11 liblapack-dev libblas-dev liblapacke-dev" >&2
  echo "or" >&2
  echo "  sudo apt install -y gfortran liblapack-dev libblas-dev liblapacke-dev" >&2
  exit 2
fi

echo "[INFO] Found all required source and mesh files."
echo "[INFO] Compiler: $FC"
"$FC" --version | head -n 1

if [[ "$FC" != *"gfortran-11"* ]]; then
  echo "[WARN] Current compiler is not gfortran-11."
  echo "[WARN] If runtime parsing differs from your local environment, try: FC=gfortran-11 ./run_cloud_compile.sh --run"
fi

echo "[INFO] Compiling buflow_test_d ..."

set -x
"$FC" -o buflow_test_d \
  TypesModule.f90 \
  meshdeformationn_d.f90 \
  BuFlow_test_d.f90 \
  main_d.f90 \
  run_parameter_d_AAA.f90 \
  -llapack
set +x

echo "[INFO] Build complete: $ROOT_DIR/buflow_test_d"

if [[ "$RUN_AFTER_BUILD" == true ]]; then
  echo "[INFO] Running: ./buflow_test_d"
  ./buflow_test_d
fi
