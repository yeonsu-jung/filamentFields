#pragma once

// Compute total entanglement by summing |lk| over provided edge pairs on GPU.
// Inputs:
//  - edges_rm: row-major flattened edges [num_edges x 6] (x1,y1,z1,x2,y2,z2)
//  - num_edges: number of edges
//  - edge_labels: label per edge (pairs with same label are skipped)
//  - pair_i, pair_j: arrays of length num_pairs specifying edge index pairs
//  - num_pairs: number of pairs
// Returns the total sum of abs(linking_number(edge_i, edge_j)).
extern "C" double ff_streaming_total_from_pairs(
    const double* edges_rm,
    int num_edges,
    const int* edge_labels,
    const int* pair_i,
    const int* pair_j,
    int num_pairs
);
