#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <math_constants.h>

extern "C"
{
__global__ void 
rtruncnorm_kernel(float *vals, int n, 
                  float *mu, float *sigma, 
                  float *lo, float *hi,
                  int rng_a, int rng_b,
                  int rng_c,
                  int maxtries)
{
    // Usual block/thread indexing...
    int myblock = blockIdx.x + blockIdx.y * gridDim.x;
    int blocksize = blockDim.x * blockDim.y * blockDim.z;
    int subthread = threadIdx.z*(blockDim.x * blockDim.y) + threadIdx.y*blockDim.x + threadIdx.x;
    int idx = myblock * blocksize + subthread;
    
    if (idx < n)
    {
            // Setup the RNG:
            curandState rng;
            curand_init (rng_a+idx*rng_b, rng_c, 0, &rng);
            // Sample:
            int accept=0;
            int numtries=0;
            while (!accept && numtries < maxtries)
            {
                  numtries ++;
                  vals[idx]=mu[idx]+sigma[idx]*curand_normal(&rng);
                  if (vals[idx]>=lo[idx] && vals[idx]<=hi[idx])
                  { accept=1; }
                  else {}
            }
                     
    }
    return;
}

} // END extern "C"