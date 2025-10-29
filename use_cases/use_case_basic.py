"""
Basic use case: create random filaments, run a single local volume analysis and print results.
"""
import numpy as np
import sys
sys.path.append('..')

import filamentFields



def make_random_filaments(n_filaments=100, points_per_filament=10, seed=0):
    rng = np.random.default_rng(seed)
    filaments = []
    for _ in range(n_filaments):
        # short random walk smoothed slightly
        x = np.cumsum(rng.standard_normal(points_per_filament))
        y = np.cumsum(rng.standard_normal(points_per_filament))
        z = np.cumsum(rng.standard_normal(points_per_filament))
        fil = np.vstack([x, y, z]).T
        filaments.append(fil)
    return filaments


def main():
    filaments = make_random_filaments(n_filaments=200, points_per_filament=12)
    fF = filamentFields.filamentFields(filaments)

    R_omega = 1.0
    query_point = np.array([0.0, 0.0, 0.0])

    # quick non-precomputed analysis
    local_edges = fF.analyze_local_volume(query_point, R_omega, 0.01)
    print("Local edges shape:", local_edges.shape)
    print("Entanglement:", fF.return_entanglement())
    print("Number of labels:", fF.return_number_of_labels())

    # precompute and run the precomputed analysis version
    fF.precompute(R_omega)
    local_edges_pc = fF.analyze_local_volume_from_precomputed(query_point, R_omega, 0.01)
    print("(precomputed) Local edges shape:", local_edges_pc.shape)
    print("(precomputed) Entanglement:", fF.return_entanglement())


if __name__ == '__main__':
    main()
