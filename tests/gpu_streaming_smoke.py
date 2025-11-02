import time
import numpy as np
import filamentFields as ff


def random_unit_vectors(n: int, rng: np.random.Generator) -> np.ndarray:
    v = rng.normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True) + 1e-12
    return v


def generate_rods(N: int, box: float, L: float, seed: int = 0):
    rng = np.random.default_rng(seed)
    centers = rng.uniform(low=-box/2, high=box/2, size=(N, 3))
    dirs = random_unit_vectors(N, rng)
    half = (L / 2.0) * dirs
    p0 = centers - half
    p1 = centers + half
    filaments = [np.vstack([p0[i], p1[i]]) for i in range(N)]
    return filaments


def main():
    N = 2000
    box = 100.0
    L = 1.0
    R = 2.0
    seed = 42

    print(f"Generating {N} rods…")
    t0 = time.time()
    filaments = generate_rods(N, box, L, seed=seed)
    t1 = time.time()
    print(f"Generation: {t1 - t0:.2f}s")

    F = ff.filamentFields(filaments)

    print("CPU streaming…")
    t2 = time.time()
    total_cpu = F.compute_total_entanglement_streaming(R)
    t3 = time.time()
    print(f"CPU: total={total_cpu:.8f}, time={t3 - t2:.3f}s")

    print("GPU streaming…")
    t4 = time.time()
    total_gpu = F.compute_total_entanglement_streaming_gpu(R)
    t5 = time.time()
    print(f"GPU: total={total_gpu:.8f}, time={t5 - t4:.3f}s")

    rel_err = abs(total_cpu - total_gpu) / max(abs(total_cpu), 1e-12)
    print(f"Relative error: {rel_err:.3e}")

    # Allow a small numerical tolerance
    ok = rel_err < 1e-9
    print("GPU streaming smoke test:", "PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
