#!/usr/bin/env bash
# SBATCH script to build with CUDA enabled and run the CUDA smoke test on FAS RC GPUs

#SBATCH -J ff_cuda_build_smoke
#SBATCH -p gpu_test
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=8000
#SBATCH --gres=gpu:1
#SBATCH -t 0-00:30
#SBATCH -o ff_cuda_build_smoke_%j.out
#SBATCH -e ff_cuda_build_smoke_%j.err

set -euo pipefail

echo "[INFO] Host: $(hostname)"
echo "[INFO] PWD:  $(pwd)"

# Load modules: let the cluster resolve default versions
module purge || true
module load cuda || true
module load cmake || true
module load gcc || true
module load python || true
echo "[INFO] Modules loaded:"
module -t list || true

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"
echo "[INFO] Repo root: $REPO_ROOT"

# Try to activate user env if available (optional)
mamba activate simdata-analysis 2>/dev/null || conda activate simdata-analysis 2>/dev/null || true

PY_EXE="$(python3 -c 'import sys; print(sys.executable)' || echo python3)"
echo "[INFO] Python: $PY_EXE"

BUILD_DIR="$REPO_ROOT/build_current"
echo "[INFO] Configure with CUDA: $BUILD_DIR"
cmake -S . -B "$BUILD_DIR" \
  -DFF_WITH_CUDA=ON \
  -DPYBIND11_FINDPYTHON=ON \
  -DPython_EXECUTABLE="$PY_EXE"

echo "[INFO] Build CUDA smoke test"
cmake --build "$BUILD_DIR" --target ff_cuda_smoke -j 4

echo "[INFO] Run CUDA smoke test"
"$REPO_ROOT/ff_cuda_smoke"

# Build the Python module as well
echo "[INFO] Build pybind11 module (filamentFields)"
cmake --build "$BUILD_DIR" --target filamentFields -j 4

# Run the Python GPU streaming smoke test
echo "[INFO] Run Python GPU streaming smoke test"
PYTHONPATH="$REPO_ROOT" "$PY_EXE" "$REPO_ROOT/tests/gpu_streaming_smoke.py"

echo "[INFO] Done"
