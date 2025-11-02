#pragma once

// Simple pair type usable from host and device
struct FFPairIJ { int x; int y; };

// Host-callable API implemented in streaming_gpu.cu
// Inputs are host pointers; the function manages device transfers internally.
extern "C" double ff_streaming_total_abs_lk_host(
    const double* h_edges_rowmajor, // length = num_edges * 6, row-major (x0,y0,z0,x1,y1,z1)
    int num_edges,
    const FFPairIJ* h_pairs,        // length = num_pairs
    int num_pairs,
    const int* h_labels             // length = num_edges
);
