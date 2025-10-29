# filamentFields — examples and use cases


This repository contains a C++ implementation with a Python pybind11 wrapper for analyzing "filament" data (collections of 3D polylines). The `use_cases/` folder contains small example scripts adapted from the original `test.py` files.

Use the examples in `use_cases/`:

 - `use_case_basic.py` — quick smoke test and single local-volume call.
 - `use_case_visualize.py` — builds a small grid of query points and saves a projected entanglement image.
 - `use_case_benchmark.py` — compares non-precomputed vs precomputed performance.

Prerequisites:

 - Python 3.8+ with `numpy` and `matplotlib`.
 - The C++ extension `filamentFields` built and discoverable by Python (see `CMakeLists.txt`).

## Build the Python extension

Use a fresh virtual environment and install basic deps (CMake, pybind11, numpy, matplotlib):

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt cmake pybind11
```

Configure and build with CMake. You can use the helper script:

```bash
PYBIND11_DIR=$(python -c 'import pybind11, pathlib;print(pathlib.Path(pybind11.__file__).resolve().parent/"share"/"cmake"/"pybind11")') \
	PYBIND11_DIR=$PYBIND11_DIR ./build_extension.sh
```

This produces `filamentFields*.so` in the repository root.

## Run examples

Run examples from the repo root with the built module on `PYTHONPATH`:

```bash
PYTHONPATH=$(pwd) python use_cases/use_case_basic.py
PYTHONPATH=$(pwd) python use_cases/use_case_visualize.py
PYTHONPATH=$(pwd) python use_cases/use_case_benchmark.py
```

`use_case_visualize.py` will save `entanglement_proj.png`.

You can also run the visualization/movie helper with a synthetic demo (no data files needed):

```bash
PYTHONPATH=$(pwd) .venv/bin/python visualize_and_compute_ent.py --demo --output-dir movie_frames_demo
```

To process your CSVs (expects files like `rod_coordinates_*.csv`), pass a positions directory:

```bash
PYTHONPATH=$(pwd) .venv/bin/python visualize_and_compute_ent.py \
	--positions-dir /path/to/positions --output-dir movie_frames --skip-frames 10 --rod-radius 1.21
```

## Run tests (smoke tests)

Install pytest and run:

```bash
pip install pytest
PYTHONPATH=$(pwd) pytest -q
```

## Notes

- If you use Anaconda, ensure CMake finds the correct pybind11 and Python headers. You can override `pybind11_DIR` via environment variable `PYBIND11_DIR`.
- The CMakeLists was updated to avoid hard-coding a specific Anaconda path for pybind11.

### Building for a specific Python version (e.g., 3.8)

Compiled extensions are specific to a Python minor version (3.8 vs 3.9 vs 3.13). To build for Python 3.8:

1) Activate a Python 3.8 environment (conda or venv)

```bash
conda create -n filfields-py38 python=3.8 -y
conda activate filfields-py38
pip install -r requirements.txt cmake pybind11
```

2) Build the extension using the active interpreter

```bash
./build_for_python.sh
```

This will produce `filamentFields.cpython-38-*.so` in the repo root. You can keep multiple versions side-by-side (e.g., 3.8 and 3.13); Python will import the one matching its version when `PYTHONPATH=$(pwd)` is used.

Note: Ensure your Python and system/toolchain architectures match (arm64 vs x86_64). On Apple Silicon, prefer arm64 Python 3.8 if available.
