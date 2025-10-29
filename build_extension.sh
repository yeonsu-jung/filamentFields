#!/usr/bin/env bash
set -euo pipefail

# Simple helper to configure and build the project with CMake.
# Run from the repository root.

BUILD_DIR=build
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
cmake ..
cmake --build . -- -j

echo "Build finished. If the Python extension was produced, run examples from the repo root or add the build dir to PYTHONPATH."
