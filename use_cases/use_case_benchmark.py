"""
Benchmark use case: compares analyze_local_volume vs analyze_local_volume_from_precomputed
over a grid of query points.
"""
import time
import numpy as np
import filamentFields


def generate_aligned_filaments(n_rods_side=30, length_points=10):
    # create a grid of short vertical rods aligned along z
    mg = np.meshgrid(np.linspace(0,1,n_rods_side), np.linspace(0,1,n_rods_side))
    mgx = mg[0].flatten()
    mgy = mg[1].flatten()
    single_edge = np.vstack([np.zeros(length_points), np.zeros(length_points), np.linspace(0,1,length_points)]).T
    filaments = []
    for i in range(len(mgx)):
        translator = np.array([mgx[i], mgy[i], 0.0])
        filaments.append(single_edge + translator)
    return filaments


def main():
    filaments = generate_aligned_filaments(n_rods_side=25, length_points=12)
    fF = filamentFields.filamentFields(filaments)

    all_nodes = np.vstack(filaments)
    xlim = (np.min(all_nodes[:,0]), np.max(all_nodes[:,0]))
    ylim = (np.min(all_nodes[:,1]), np.max(all_nodes[:,1]))
    zlim = (np.min(all_nodes[:,2]), np.max(all_nodes[:,2]))

    num_grid = 15
    mg = np.meshgrid(np.linspace(xlim[0], xlim[1], num_grid),
                     np.linspace(ylim[0], ylim[1], num_grid),
                     np.linspace(zlim[0], zlim[1], num_grid))
    query_points = np.vstack([mg[0].flatten(), mg[1].flatten(), mg[2].flatten()]).T

    R_omega = 0.05

    # benchmark non-precomputed
    t0 = time.time()
    ent_non = np.zeros(len(query_points))
    for i, qp in enumerate(query_points):
        fF.analyze_local_volume(qp, R_omega, 0.01)
        ent_non[i] = fF.return_entanglement()
    dt_non = time.time() - t0

    # benchmark precomputed
    fF.precompute(R_omega)
    t1 = time.time()
    ent_pc = np.zeros(len(query_points))
    for i, qp in enumerate(query_points):
        fF.analyze_local_volume_from_precomputed(qp, R_omega, 0.01)
        ent_pc[i] = fF.return_entanglement()
    dt_pc = time.time() - t1

    print(f'Non-precomputed time: {dt_non:.3f} s')
    print(f'Precomputed time:     {dt_pc:.3f} s')
    print(f'Average entanglement (non): {np.nanmean(ent_non):.6f}')
    print(f'Average entanglement (pc):  {np.nanmean(ent_pc):.6f}')


if __name__ == '__main__':
    main()
