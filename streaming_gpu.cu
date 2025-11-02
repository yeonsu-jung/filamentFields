#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <cmath>
#include "cuda_streaming_api.cuh"

#define INV_TWO_PI 0.1591549430918953357688837633725143620344596457404564487476673440

__device__ inline void load_edge_rm(const double* edges_rm, int idx, double e[6]) {
    // Row-major layout: 6 consecutive doubles per edge
    const double* p = edges_rm + idx * 6;
    #pragma unroll
    for (int k = 0; k < 6; ++k) e[k] = p[k];
}

__device__ inline void vec3_sub(const double a[3], const double b[3], double out[3]) {
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
}

__device__ inline double vdot(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

__device__ inline void vcross(const double a[3], const double b[3], double out[3]) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

__device__ inline double vnorm(const double a[3]) {
    return sqrt(vdot(a,a));
}

__device__ inline double compute_lk_for_edges_rm(const double e1[6], const double e2[6]) {
    // Match CPU formula using Qu and James variant
    const double p0[3] = {e1[0], e1[1], e1[2]};
    const double p1[3] = {e1[3], e1[4], e1[5]};
    const double q0[3] = {e2[0], e2[1], e2[2]};
    const double q1[3] = {e2[3], e2[4], e2[5]};

    double a[3], b[3], c[3], d[3];
    // a = p0 - q0; b = p0 - q1; c = p1 - q1; d = p1 - q0
    vec3_sub(p0, q0, a);
    vec3_sub(p0, q1, b);
    vec3_sub(p1, q1, c);
    vec3_sub(p1, q0, d);

    double bxc[3], dxa[3];
    vcross(b, c, bxc);
    vcross(d, a, dxa);

    double na = vnorm(a), nb = vnorm(b), nc = vnorm(c), nd = vnorm(d);

    double num1 = vdot(a, bxc);
    double den1 = na*nb*nc + vdot(a,b)*nc + vdot(c,a)*nb + vdot(b,c)*na;

    double num2 = vdot(c, dxa);
    double den2 = nc*nd*na + vdot(c,d)*na + vdot(a,c)*nd + vdot(d,a)*nc;

    double lk = INV_TWO_PI * (atan2(num1, den1) + atan2(num2, den2));
    return lk;
}

__global__ void kernel_pairs_abs_lk(
    const double* edges_rm, int num_edges,
    const FFPairIJ* pairs, int num_pairs,
    const int* labels,
    double* out_abs
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_pairs) return;
    int idx = pairs[i].x;
    int jdx = pairs[i].y;
    if (idx == jdx) { out_abs[i] = 0.0; return; }
    if (labels && labels[idx] == labels[jdx]) { out_abs[i] = 0.0; return; }
    double e1[6], e2[6];
    load_edge_rm(edges_rm, idx, e1);
    load_edge_rm(edges_rm, jdx, e2);
    double lk = compute_lk_for_edges_rm(e1, e2);
    out_abs[i] = fabs(lk);
}

extern "C" double ff_streaming_total_abs_lk_host(
    const double* h_edges_rowmajor, int num_edges,
    const FFPairIJ* h_pairs, int num_pairs,
    const int* h_labels
) {
    if (num_pairs <= 0 || num_edges <= 0) return 0.0;

    // Device buffers
    double *d_edges = nullptr, *d_abs = nullptr;
    FFPairIJ* d_pairs = nullptr;
    int* d_labels = nullptr;

    size_t edges_bytes = static_cast<size_t>(num_edges) * 6 * sizeof(double);
    size_t pairs_bytes = static_cast<size_t>(num_pairs) * sizeof(FFPairIJ);
    cudaMalloc(&d_edges, edges_bytes);
    cudaMalloc(&d_pairs, pairs_bytes);
    cudaMalloc(&d_abs, static_cast<size_t>(num_pairs) * sizeof(double));
    if (h_labels) cudaMalloc(&d_labels, static_cast<size_t>(num_edges) * sizeof(int));

    cudaMemcpy(d_edges, h_edges_rowmajor, edges_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_pairs, h_pairs, pairs_bytes, cudaMemcpyHostToDevice);
    if (h_labels) cudaMemcpy(d_labels, h_labels, static_cast<size_t>(num_edges) * sizeof(int), cudaMemcpyHostToDevice);

    int block = 256;
    int grid = (num_pairs + block - 1) / block;
    kernel_pairs_abs_lk<<<grid, block>>>(d_edges, num_edges, d_pairs, num_pairs, d_labels, d_abs);
    cudaError_t err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        // Cleanup and throw value indicating failure (caller can detect negative)
        cudaFree(d_edges); cudaFree(d_pairs); cudaFree(d_abs); if (d_labels) cudaFree(d_labels);
        return -1.0;
    }

    // Reduce on device
    thrust::device_ptr<double> dptr(d_abs);
    double total = thrust::reduce(dptr, dptr + num_pairs, 0.0, thrust::plus<double>());

    cudaFree(d_edges);
    cudaFree(d_pairs);
    cudaFree(d_abs);
    if (d_labels) cudaFree(d_labels);
    return total;
}
