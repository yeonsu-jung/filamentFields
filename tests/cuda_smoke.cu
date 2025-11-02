#include <cstdio>
#include <cuda_runtime.h>

__global__ void add_one(const int* in, int* out, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = in[i] + 1;
}

int main() {
    const int N = 256;
    int h_in[N], h_out[N];
    for (int i = 0; i < N; ++i) h_in[i] = i;

    int *d_in = nullptr, *d_out = nullptr;
    cudaError_t err;

    err = cudaMalloc(&d_in, N * sizeof(int));
    if (err != cudaSuccess) { std::fprintf(stderr, "cudaMalloc d_in failed: %s\n", cudaGetErrorString(err)); return 1; }
    err = cudaMalloc(&d_out, N * sizeof(int));
    if (err != cudaSuccess) { std::fprintf(stderr, "cudaMalloc d_out failed: %s\n", cudaGetErrorString(err)); return 1; }

    err = cudaMemcpy(d_in, h_in, N * sizeof(int), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) { std::fprintf(stderr, "cudaMemcpy H->D failed: %s\n", cudaGetErrorString(err)); return 1; }

    dim3 block(128);
    dim3 grid((N + block.x - 1) / block.x);
    add_one<<<grid, block>>>(d_in, d_out, N);
    err = cudaDeviceSynchronize();
    if (err != cudaSuccess) { std::fprintf(stderr, "Kernel failed: %s\n", cudaGetErrorString(err)); return 1; }

    err = cudaMemcpy(h_out, d_out, N * sizeof(int), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) { std::fprintf(stderr, "cudaMemcpy D->H failed: %s\n", cudaGetErrorString(err)); return 1; }

    // Verify
    for (int i = 0; i < 5; ++i) {
        std::printf("%d -> %d\n", h_in[i], h_out[i]);
    }
    bool ok = true;
    for (int i = 0; i < N; ++i) if (h_out[i] != h_in[i] + 1) { ok = false; break; }
    std::printf("CUDA smoke test: %s\n", ok ? "PASS" : "FAIL");

    cudaFree(d_in);
    cudaFree(d_out);
    return ok ? 0 : 1;
}
