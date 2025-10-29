"""
Visualize a set of filaments and compute an entanglement field on a small grid.
"""
import numpy as np
import matplotlib.pyplot as plt
import filamentFields


def make_smoothed_random_filaments(n_filaments=200, points=100, seed=1):
    rng = np.random.default_rng(seed)
    filaments = []
    for _ in range(n_filaments):
        x = np.cumsum(rng.standard_normal(points))
        y = np.cumsum(rng.standard_normal(points))
        z = np.cumsum(rng.standard_normal(points))
        # simple moving average smoothing
        k = 5
        x = np.convolve(x, np.ones(k)/k, mode='valid')
        y = np.convolve(y, np.ones(k)/k, mode='valid')
        z = np.convolve(z, np.ones(k)/k, mode='valid')
        filaments.append(np.vstack([x, y, z]).T)
    return filaments


def main():
    filaments = make_smoothed_random_filaments(n_filaments=300, points=60)
    fF = filamentFields.filamentFields(filaments)

    all_nodes = np.vstack(filaments)
    xlim = (np.min(all_nodes[:,0]), np.max(all_nodes[:,0]))
    ylim = (np.min(all_nodes[:,1]), np.max(all_nodes[:,1]))
    zlim = (np.min(all_nodes[:,2]), np.max(all_nodes[:,2]))

    # small grid for visualization
    num_grid = 10
    mg = np.meshgrid(np.linspace(xlim[0], xlim[1], num_grid),
                     np.linspace(ylim[0], ylim[1], num_grid),
                     np.linspace(zlim[0], zlim[1], num_grid))
    query_points = np.vstack([mg[0].flatten(), mg[1].flatten(), mg[2].flatten()]).T

    R_omega = max(np.ptp(all_nodes[:,0]), np.ptp(all_nodes[:,1]), np.ptp(all_nodes[:,2])) * 0.1

    # use precompute for speed
    fF.precompute(R_omega)
    results = fF.analyze_local_volume_over_domain_from_precomputed(query_points, R_omega, 0.01)

    entanglement = results[:,3].reshape(num_grid, num_grid, num_grid)
    # project along the z axis for a 2D image
    ent_proj = np.sum(entanglement, axis=2)

    plt.figure()
    plt.imshow(ent_proj)
    plt.title('Projected entanglement (z)')
    plt.colorbar()
    plt.savefig('entanglement_proj.png')
    print('Saved entanglement_proj.png')


if __name__ == '__main__':
    main()
