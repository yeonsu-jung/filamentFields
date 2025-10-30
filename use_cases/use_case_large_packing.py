"""
Large packing benchmark: generate N random rods (2-node open segments) with random
positions and orientations, then compute global entanglement via streaming (no
NxN matrix allocation).

Usage:
  PYTHONPATH=. python use_cases/use_case_large_packing.py --N 100000 --box 100.0 --L 1.0 --R 2.0

Tips:
- Choose R (interaction radius) relative to rod length L; larger R means more candidate pairs.
- 100k rods may take time depending on R and box size; try N=20000 first to gauge speed.
"""
import argparse
import time
import numpy as np
import filamentFields


def random_unit_vectors(n: int, rng: np.random.Generator) -> np.ndarray:
    v = rng.normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True) + 1e-12
    return v


def generate_rods(N: int, box: float, L: float, seed: int = 0):
    """Return list of N rods, each as a (2,3) array of endpoints.
    Centers are uniform in [-box/2, box/2]^3, orientations random, length L.
    """
    rng = np.random.default_rng(seed)
    centers = rng.uniform(low=-box/2, high=box/2, size=(N, 3))
    dirs = random_unit_vectors(N, rng)
    half = (L / 2.0) * dirs
    p0 = centers - half
    p1 = centers + half
    filaments = [np.vstack([p0[i], p1[i]]) for i in range(N)]
    return filaments


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=100_000, help='Number of rods (filaments)')
    ap.add_argument('--box', type=float, default=100.0, help='Box size (cube side length)')
    ap.add_argument('--L', type=float, default=1.0, help='Rod length')
    ap.add_argument('--R', type=float, default=2.0, help='Interaction radius R_omega for pairing')
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()

    print(f"Generating {args.N:,} rods (L={args.L}, box={args.box})…")
    t0 = time.time()
    filaments = generate_rods(args.N, args.box, args.L, seed=args.seed)
    t1 = time.time()
    print(f"Generation: {(t1 - t0):.2f}s")

    print("Building filamentFields object…")
    fF = filamentFields.filamentFields(filaments)

    print(f"Computing streaming total entanglement with R={args.R}…")
    t2 = time.time()
    total = fF.compute_total_entanglement_streaming(args.R)
    t3 = time.time()

    pairs = fF.return_edge_pair_count()
    print(f"Candidate pairs: {pairs:,}")
    print(f"Total entanglement (sum |lk|): {total:.6f}")
    print(f"Compute time: {(t3 - t2):.2f}s (excluding generation)")


if __name__ == '__main__':
    main()
