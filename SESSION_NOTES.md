# Session Notes: filamentFields acceleration and large-scale use

This file summarizes the key changes and how to reproduce results on another machine (e.g., a cluster).

## Features implemented

- Broad-phase AABB sweep for edge pairing (CPU, no GPU required)
- Streaming global entanglement: exact Gauss segment kernel without NxN matrix
  - Python: `compute_total_entanglement_streaming(R_omega)`
- Barnes–Hut (BH) near/far split (CPU) with opening angle `theta`
  - Python: `compute_total_entanglement_bh(theta=0.5, maxLeaf=64)` returns `(total, error_bound)`
- Large packing use cases
  - `use_cases/use_case_large_packing.py`: synthetic N rods
  - `use_cases/use_case_large_csv.py`: load `i,cx,cy,cz,phi,theta` CSV (e.g., `tests/cfg_placed_121600.csv`)

## How to set up on a new machine (cluster)

1. Clone the repo

```zsh
git clone https://github.com/yeonsu-jung/filamentFields.git
cd filamentFields
```

1. Create/activate a Python environment (example with conda, Python 3.12)

```zsh
conda create -n ff-env python=3.12 -y
conda activate ff-env
pip install -r requirements.txt
pip install pybind11 cmake
```

1. Build the extension for the active Python

Option A: helper script

```zsh
bash ./build_for_python.sh
```

Option B: manual cmake

```zsh
cmake -S . -B build_current \
  -DPYBIND11_FINDPYTHON=ON \
  -DPython_EXECUTABLE="$(python -c 'import sys; print(sys.executable)')"
cmake --build build_current -j 8
```

1. Sanity check


```zsh
python - <<'PY'
import filamentFields as ff
print('OK:', ff)
PY
```

## CUDA/GPU build and test on FAS RC

On the FAS RC cluster, compile CUDA code on a GPU node and load the CUDA toolkit module.

- Interactive session (quick):

```bash
salloc -p gpu_test -t 0-01:00 --mem 8000 --gres=gpu:1
module load cuda
module load cmake
# Optional: module load gcc
cd /path/to/rod-placement/extern/filamentFields
cmake -S . -B build_current -DFF_WITH_CUDA=ON \
  -DPYBIND11_FINDPYTHON=ON \
  -DPython_EXECUTABLE="$(python3 -c 'import sys; print(sys.executable)')"
cmake --build build_current --target ff_cuda_smoke -j 8
./ff_cuda_smoke
```

- Batch job (recommended repeatable):

```bash
cd /path/to/rod-placement/extern/filamentFields
sbatch sbatch_cuda_build_and_smoke.sh
```

Notes:
- SLURM GPUs: request with `--gres=gpu:1` and pick a GPU partition such as `gpu_test`, `gpu`, or `gpu_requeue`.
- You can check the device and driver/toolkit with `nvidia-smi` on the GPU node.
- Default CUDA architectures in CMake are set to build for V100/A100/H100 (70;80;90). Override with `-DCMAKE_CUDA_ARCHITECTURES=80` if desired.

### Build the Python module with CUDA and run the GPU streaming test

```bash
# on a GPU node with cuda + cmake modules loaded
cd /path/to/rod-placement/extern/filamentFields
cmake -S . -B build_current -DFF_WITH_CUDA=ON \
  -DPYBIND11_FINDPYTHON=ON \
  -DPython_EXECUTABLE="$(python3 -c 'import sys; print(sys.executable)')"
cmake --build build_current -j 8

# Run the Python GPU/CPU comparison
PYTHONPATH=. python tests/gpu_streaming_smoke.py
```

This compares CPU `compute_total_entanglement_streaming(R)` with GPU `compute_total_entanglement_streaming_gpu(R)` on a small synthetic set and prints timings and relative error.


## Running the large CSV (121,600 rods)

- Streaming (R-limited) global total:

```zsh
PYTHONPATH=. python use_cases/use_case_large_csv.py \
  --csv tests/cfg_placed_121600.csv --L 0.2 --R 0.3
```

- BH (full total with error bound):

```zsh
python - <<'PY'
import numpy as np
import filamentFields as ff
from use_cases.use_case_large_csv import load_cfg_csv, rods_from_centers_dirs
centers, dirs = load_cfg_csv('tests/cfg_placed_121600.csv')
filaments = rods_from_centers_dirs(centers, dirs, L=0.2)
F = ff.filamentFields(filaments)
# Tune theta (0.6–0.8 for speed, 0.4–0.6 for accuracy) and maxLeaf 64–256
total, err = F.compute_total_entanglement_bh(theta=0.7, maxLeaf=128)
print('BH total=', total, 'err=', err, 'rel_err~', err/max(total,1e-12))
PY
```

## Notes

- For exact (all-pairs) totals at small N, `compute_total_entanglement_bh(theta=0.5)` will typically refine most nodes, yielding zero error bound.
- For R-limited totals, prefer the streaming API (BH targets near/far accuracy rather than a geometric cutoff).
- On macOS, GPU acceleration isn’t available via CUDA; plan to run GPU kernels on a Linux/NVIDIA cluster later.

## What changed in code

- `filamentFields.h/.cpp`: AABB pairing, streaming method, BH octree/build/traverse, and error bound accumulation.
- `wrapper.cpp`: exposed the new methods to Python.
- `use_cases`: added `use_case_large_packing.py` and `use_case_large_csv.py`.

## Recreating our quick benchmarks

- Synthetic 20k rods, L=1, box=100, R=2 (streaming):

```zsh
PYTHONPATH=. python use_cases/use_case_large_packing.py --N 20000 --box 100 --L 1.0 --R 2.0
```

- 121,600 rods CSV (streaming with L=0.2, R=0.3): ~216M candidate pairs, ~20–30s on a modern CPU.

## Optional: archive this repo with history (offline transfer)

```zsh
git bundle create filamentFields.bundle --all
# move the .bundle file to the cluster
# then on cluster:
git clone filamentFields.bundle filamentFields
```

## Optional: capture exact Python deps

```zsh
pip freeze > pinned-requirements.txt
```

This document should be enough to resume work on a fresh machine and reproduce the results and methods added in this session.
