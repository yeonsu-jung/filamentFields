"""
Compute streaming total entanglement for a large packing given as a CSV with columns:
  i,cx,cy,cz,phi,theta
Each row defines a rod centered at (cx,cy,cz) with orientation given by spherical
angles (phi, theta). Endpoints are computed for a chosen rod length L.

Usage:
  PYTHONPATH=. python use_cases/use_case_large_csv.py \
      --csv tests/cfg_placed_121600.csv --L 0.2 --R 0.3

Notes:
- Orientation mapping uses spherical convention: dir = [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)].
- If your convention differs, tweak the mapping below accordingly.
- For very large R relative to point spacing and L, candidate pairs can grow quickly.
"""
import argparse
import time
import numpy as np
import filamentFields


def load_cfg_csv(path: str):
    data = np.loadtxt(path, delimiter=',', skiprows=1)
    # columns: idx, cx, cy, cz, phi, theta
    centers = data[:, 1:4]
    phi = data[:, 4]
    theta = data[:, 5]
    # spherical to Cartesian unit vector
    s = np.sin(theta)
    c = np.cos(theta)
    dirs = np.stack([np.cos(phi) * s, np.sin(phi) * s, c], axis=1)
    # normalize just in case
    dirs /= (np.linalg.norm(dirs, axis=1, keepdims=True) + 1e-12)
    return centers, dirs


def rods_from_centers_dirs(centers: np.ndarray, dirs: np.ndarray, L: float):
    half = (L / 2.0) * dirs
    p0 = centers - half
    p1 = centers + half
    filaments = [np.vstack([p0[i], p1[i]]) for i in range(centers.shape[0])]
    return filaments


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True, help='Path to cfg_placed_*.csv with columns i,cx,cy,cz,phi,theta')
    ap.add_argument('--L', type=float, default=0.2, help='Rod length')
    ap.add_argument('--R', type=float, default=0.3, help='Interaction radius R_omega for pairing')
    args = ap.parse_args()

    t0 = time.time()
    centers, dirs = load_cfg_csv(args.csv)
    t1 = time.time()
    print(f"Loaded {centers.shape[0]:,} rods from {args.csv} in {(t1-t0):.2f}s")

    filaments = rods_from_centers_dirs(centers, dirs, args.L)
    t2 = time.time()
    print(f"Constructed filament list in {(t2-t1):.2f}s")

    fF = filamentFields.filamentFields(filaments)
    print(f"Computing streaming total entanglement with R={args.R}…")
    t3 = time.time()
    total = fF.compute_total_entanglement_streaming(args.R)
    t4 = time.time()

    pairs = fF.return_edge_pair_count()
    print(f"Candidate pairs: {pairs:,}")
    print(f"Total entanglement (sum |lk|): {total:.6f}")
    print(f"Compute time: {(t4 - t3):.2f}s (excluding IO and construction)")


if __name__ == '__main__':
    main()
