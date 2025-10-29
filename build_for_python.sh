#!/usr/bin/env bash
set -euo pipefail

# Build the pybind11 extension for the CURRENT Python interpreter.
# Usage:
#   1) Activate the Python you want (e.g., conda/env with Python 3.8)
#   2) pip install -r requirements.txt cmake pybind11
#   3) ./build_for_python.sh

PY_EXE=${PY_EXE:-$(python -c 'import sys; print(sys.executable)')}
echo "Using Python: $PY_EXE"

# Optional: help pybind11 find the right CMake package
if [[ -z "${PYBIND11_DIR:-}" ]]; then
  PYBIND11_DIR=$("$PY_EXE" - <<'PY'
import pybind11, pathlib
print(pathlib.Path(pybind11.__file__).resolve().parent/"share"/"cmake"/"pybind11")
PY
  )
fi
echo "Using PYBIND11_DIR: $PYBIND11_DIR"

BUILD_DIR=build_current
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake .. \
  -DPYBIND11_FINDPYTHON=ON \
  -DPython_EXECUTABLE="$PY_EXE" \
  -DPybind11_DIR="$PYBIND11_DIR"

cmake --build . -- -j

echo "Done. If successful, a filamentFields*.so matching this Python will be in the repository root."
