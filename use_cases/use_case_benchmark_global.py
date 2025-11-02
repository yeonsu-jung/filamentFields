"""
Global entanglement benchmark: compare performance and accuracy of
- streaming (R-limited or full with very large R)
- BH near/far with opening angle theta

Usage examples:
  # Moderate N with baseline = full streaming (very large R)
  PYTHONPATH=. python use_cases/use_case_benchmark_global.py \
      --N 5000 --box 100 --L 1.0 --baseline-R 1e9 \
      --thetas 0.4 0.6 0.8 --maxLeaf 128

  # Large N (skip baseline): compare BH vs streaming with finite R
  PYTHONPATH=. python use_cases/use_case_benchmark_global.py \
      --N 50000 --box 100 --L 1.0 --stream-R 2.0 \
      --thetas 0.6 0.8 --maxLeaf 128 --skip-baseline
"""
from __future__ import annotations
import argparse
import time
import numpy as np
import filamentFields as ff
from typing import Iterable
from use_cases.use_case_large_csv import load_cfg_csv, rods_from_centers_dirs


def random_unit_vectors(n: int, rng: np.random.Generator) -> np.ndarray:
    v = rng.normal(size=(n, 3))
    v /= (np.linalg.norm(v, axis=1, keepdims=True) + 1e-12)
    return v


def generate_rods(N: int, box: float, L: float, seed: int = 0):
    rng = np.random.default_rng(seed)
    centers = rng.uniform(-box/2, box/2, size=(N, 3))
    dirs = random_unit_vectors(N, rng)
    half = (L / 2.0) * dirs
    p0 = centers - half
    p1 = centers + half
    filaments = [np.vstack([p0[i], p1[i]]) for i in range(N)]
    return filaments


def time_call(fn, *args, **kwargs):
    t0 = time.perf_counter()
    out = fn(*args, **kwargs)
    t1 = time.perf_counter()
    return out, t1 - t0


def fmt_s(x: float) -> str:
    return f"{x:.3f}s"


def main():
    ap = argparse.ArgumentParser()
    # Data source: synthetic or CSV
    ap.add_argument('--csv', type=str, default=None, help='Path to cfg_placed_*.csv (i,cx,cy,cz,phi,theta). If set, CSV is used instead of synthetic.')
    ap.add_argument('--N', type=int, default=5000, help='Number of rods for synthetic generation')
    ap.add_argument('--box', type=float, default=100.0, help='Box size for synthetic generation')
    ap.add_argument('--L', type=float, default=1.0, help='Rod length (used for both synthetic and CSV)')
    ap.add_argument('--seed', type=int, default=0)

    # Streaming params
    ap.add_argument('--baseline-R', type=float, default=1e9, help='R for baseline streaming (very large => all pairs)')
    ap.add_argument('--stream-R', type=float, default=2.0, help='R for additional streaming test (optional)')
    ap.add_argument('--skip-baseline', action='store_true', help='Skip full streaming baseline (useful for very large N)')

    # BH params
    ap.add_argument('--thetas', type=float, nargs='*', default=[0.6, 0.8], help='Opening angles to test')
    ap.add_argument('--maxLeaf', type=int, default=128, help='Max leaf size for BH tree')

    args = ap.parse_args()

    if args.csv:
        print(f"Loading rods from CSV: {args.csv} (L={args.L})…")
        centers, dirs = load_cfg_csv(args.csv)
        filaments = rods_from_centers_dirs(centers, dirs, args.L)
        print(f"Loaded {len(filaments):,} rods from CSV.")
    else:
        print(f"Generating N={args.N:,} rods, L={args.L}, box={args.box} (seed={args.seed})…")
        filaments = generate_rods(args.N, args.box, args.L, seed=args.seed)
    F = ff.filamentFields(filaments)

    # Baseline: streaming with very large R (effectively all pairs)
    baseline_total = None
    if not args.skip_baseline:
        print(f"Baseline (streaming) with R={args.baseline_R}…")
        (baseline_total,), t_base = time_call(lambda: (F.compute_total_entanglement_streaming(args.baseline_R),))
        print(f"  total = {baseline_total:.6f} | time = {fmt_s(t_base)}")

    # Optional: streaming with finite R for comparison
    if args.stream_R is not None and args.stream_R > 0:
        print(f"Streaming with R={args.stream_R}…")
        (stream_total,), t_stream = time_call(lambda: (F.compute_total_entanglement_streaming(args.stream_R),))
        if baseline_total is not None:
            diff = abs(stream_total - baseline_total)
            rel = diff / max(abs(baseline_total), 1e-12)
            print(f"  total = {stream_total:.6f} | time = {fmt_s(t_stream)} | Δ={diff:.6f} (rel {100*rel:.3f}%)")
        else:
            print(f"  total = {stream_total:.6f} | time = {fmt_s(t_stream)}")

    # BH sweeps
    for theta in args.thetas:
        print(f"BH(theta={theta}, maxLeaf={args.maxLeaf})…")
        (bh_total, bh_err), t_bh = time_call(F.compute_total_entanglement_bh, theta, args.maxLeaf)
        if baseline_total is not None:
            diff = abs(bh_total - baseline_total)
            rel = diff / max(abs(baseline_total), 1e-12)
            print(f"  total = {bh_total:.6f} | bound = {bh_err:.6f} | time = {fmt_s(t_bh)} | Δ={diff:.6f} (rel {100*rel:.3f}%)")
        else:
            print(f"  total = {bh_total:.6f} | bound = {bh_err:.6f} | time = {fmt_s(t_bh)}")


if __name__ == '__main__':
    main()
