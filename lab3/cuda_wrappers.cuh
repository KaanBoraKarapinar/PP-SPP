#pragma once
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <iostream>

__host__ cudaError_t parprog_cudaMemcpy( void* dst, const void* src, size_t count, cudaMemcpyKind kind ){
   return cudaMemcpy(dst, src, count, kind);
}

__host__ cudaError_t parprog_cudaMalloc( void** devPtr, size_t size ){
   return cudaMalloc(devPtr, size);
}

__host__  cudaError_t parprog_cudaFree( void* devPtr ){
   return cudaFree(devPtr);
}