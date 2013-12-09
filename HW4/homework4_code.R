##Question 1
##a
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
    int subthread = threadIdx.z*(blockDim.x * blockDim.y) 
    + threadIdx.y*blockDim.x + threadIdx.x;
    
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
      if (numtries==maxtries)
      {
        // Handle the One Sided Truncated Normal, 
        Distribution Tail Case:
          if (isfinite(hi[idx])==0)
          {
            float u_lo=(lo[idx]-mu[idx])/sigma[idx];
            float alpha=(u_lo+sqrt((u_lo*u_lo)+4))/2;
            int accept2=0;
            while (!accept2)
            {
              float z=u_lo-log(curand_uniform(&rng))/alpha;
              float u1=curand_uniform(&rng);
              float qz;
              if (u_lo < alpha)
              {qz=exp(-((alpha-z)*(alpha-z))/2);}
              else {qz=exp(-((u_lo-alpha)*
                               (u_lo-alpha))/2-((alpha-z)*(alpha-z))/2);}
              if (u1<qz) 
              {
                accept2=1;
                vals[idx]=sigma[idx]*z+mu[idx];
              }
              else {}
            }
          }
        
        else if  (isfinite(lo[idx])==0)
        {
          float u_lo=-(hi[idx]-mu[idx])/sigma[idx];
          float alpha=(u_lo+sqrt((u_lo*u_lo)+4))/2;
          int accept3=0;
          while (!accept3)
          {
            float z=u_lo
            -log(curand_uniform(&rng))/alpha;
            float u1=curand_uniform(&rng);
            float qz;
            if (u_lo < alpha)
            {qz=exp(-((alpha-z)*(alpha-z))/2);}
            else {qz=exp(-((u_lo-alpha)*(u_lo-alpha))
                         /2-((alpha-z)*(alpha-z))/2);}
            if (u1<qz) 
            {
              accept3=1;
              vals[idx]=-sigma[idx]*z+mu[idx];
            }
            else {}        
          }
        }
        // Handle the Two Sided Truncated Normal, 
        Distribution Tail Case:
          else 
          {
            int accept4=0;
            while (! accept4)
            {
              float z;
              float q;
              float u;
              z=lo[idx]+(hi[idx]-lo[idx])*
                curand_uniform(&rng);
              
              if (0>=lo[idx] && 0<=hi[idx])
              {
                q=exp(-(z*z)/2) ;
              }
              else if (hi[idx]<0)
              {
                q=exp((hi[idx]*hi[idx]-z*z)/2);
              }
              else
              {
                q=exp((lo[idx]*lo[idx]-z*z)/2);
              }
              u=curand_uniform(&rng);
              if (u<=q)
              {
                accept4=1;
                vals[idx]=z;
              }
              else {}
            }
          }                 
      }                     
    }
    return;
  }
} // END extern "C"

##b
rcuda_trun=function(N, rng_a, rng_b, rng_c, maxtries)
{
  library(RCUDA)
  cat("Loading module...\n") 
  m = loadModule("truncnorm.ptx")
  
  cat("done. Extracting kernels...\n")
  k_rtruncnorm = m$rtruncnorm_kernel
  
  cat("done. Setting up Parameters...\n")
  ##Maybe revised when mu, sigma, lo, hi are not the same
  lo = rep(0, N)
  hi = rep(1.5, N)
  mu =rep(2, N)
  sigma=rep(1, N)
  
  cat("done. Calculating Grid/ Block dimensions...\n")
  bg = compute_grid(N)
  grid_dims = bg$grid_dims
  block_dims = bg$block_dims
  
  cat("Grid size:\n")
  print(grid_dims)
  cat("Block size:\n")
  print(block_dims)
  
  nthreads = prod(grid_dims)*prod(block_dims) 
  cat("Total number of threads to launch = ",nthreads,"\n")
  
  cat("done. Runing CUDA kernels...\n")
  x = rep(0,N)
  all_time=system.time({
    mem = copyToDevice(x) 
    .cuda(k_rtruncnorm, mem, N, mu, sigma, lo, hi, 
          rng_a, rng_b, rng_c, maxtries, 
          gridDim=grid_dims, blockDim=block_dims)
    
    cat("Copying result back from device...\n")
    trun_num = copyFromDevice(obj=mem, nels=mem$nels, type="float")
  })
  
  write.table(trun_num, "trun.txt", sep=" ")
  cat("Time: \n")
  print(all_time)
  cat("Mean: \n")
  print(mean(trun_num))
  return(all_time)
}

##d
#Inverse CDF
trun_2=function (mu, sigma, a, b, n)
{
  num=runif(n, min=0, max=1)
  alpha=(a-mu)/sigma
  beta=(b-mu)/sigma
  sam=sapply(num, function(x) qnorm(pnorm(alpha)+x*(pnorm(beta)
                                                    -pnorm(alpha)))*sigma+mu)
  return(sam)
}

#Alternative version (Same algorithm with the GPU one)
trun = function(mu, sigma, a, b, maxtries=2000)
{
  accept=FALSE
  numtries=0
  
  while(!accept && numtries < maxtries)
  {
    numtries=numtries+1
    x=rnorm(1, mu, sigma)
    if (x>=a && x<=b)
    {
      accept=TRUE
    }
    else {}
  }
  
  if (numtries==maxtries)
  {
    if (b==Inf)
    {
      u_lo=(a-mu)/sigma
      alpha=(u_lo+sqrt((u_lo^2)+4))/2
      accept2=FALSE
      while (!accept2)
      {
        z=u_lo+rexp(1, alpha)
        u=runif(1, 0, 1)
        if (u_lo < alpha) {qz=exp(-((alpha-z)^2)/2)}
        else {qz=exp(-((u_lo-alpha)^2)/2-((alpha-z)^2)/2)}
        if (u<qz) 
        {accept2=TRUE
         x=sigma*z+mu
        }
        else
        {}
      }
    }
    
    if (a==-Inf)
    {
      u_lo=-(b-mu)/sigma
      alpha=(u_lo+sqrt((u_lo^2)+4))/2
      accept3=FALSE
      while (!accept3)
      {
        z=u_lo+rexp(1, alpha)
        u=runif(1, 0, 1)
        if (u_lo < alpha) {qz=exp(-((alpha-z)^2)/2)}
        else {qz=exp(-((u_lo-alpha)^2)/2-((alpha-z)^2)/2)}
        if (u<qz) 
        {accept3=TRUE
         x=-sigma*z+mu
        }
        else
        {}
      }
    }
  }
  return(x)
}

##e
n_vec=as.integer(c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8))

#System time for RCUDA function:
time2=sapply(n_vec, function(n) 
  rcuda_trun(n, rng_a, rng_b, rng_c, maxtries))

#System time for pure R code:
#Need to repeat a couple of times if the n is to small to calculate 
#running time
time=sapply(n_vec,
            function(n){
              if(n<100)
                B=10000
              else if(n<1000)
                B=100
              else if (n<10000)
                B=10
              else
                B=1
              system.time(replicate(B,trun_2(2, 1, 0, 1.5, n)))/B
            })

##Question 2
##a
probit_mcmc_cpu = function(y, X, 
                           p,  ##length of beta
                           N,  ##Number of Observations
                           sigma,##sigma for truncated normal z
                           lo, ##lower bound for truncated normal z
                           hi, ##upper bound for truncated normal z 
                           beta.0, ##beta_prior mean
                           Sigma.0.inv, ## beta prior variance inverse
                           niter=2000, burnin=500)
{
  beta=matrix(numeric(p*(burnin+niter+1)), (burnin+niter+1), p) 
  for (i in 1: (burnin+niter))
  {
    ##Sample from truncated normal for z
    u=as.vector(X%*%beta[i, ])
    z=rtruncnorm(1, lo, hi, u, sigma)
    
    ##Calculate beta_posterior mean and variance
    po_beta_mu=solve(Sigma.0.inv+t(X)%*%X)%*%
      (Sigma.0.inv%*%beta.0+t(X)%*%z)
    po_beta_sigma=solve(Sigma.0.inv+t(X)%*%X)
    
    ##Sample beta from multivariate normal
    beta[i+1,]=mvrnorm(n = 1, mu=po_beta_mu, Sigma=po_beta_sigma)
    
    ##Print results for every 250 iterations
    if(i%%250==0)
    {
      cat("Finished", i, "Iterations...\n")
    }
  }
  cat("Work Completed...\n")
  return(beta[(burnin+2):(burnin+niter+1),])
}

##b
probit_mcmc_gpu=function(y, X, 
                         p, ##length of beta
                         N, ##Number of Observations
                         sigma, ##sigma for truncated normal z
                         val, ##parameter to store the simulated z
                         lo, ##lower bound for truncated normal z
                         hi, ##upper bound for truncated normal z
                         rng_a, rng_b, ##RNG parameter
                         maxtries, 
                         m, k_rtruncnorm, 
                         beta_0, ##beta_prior mean
                         Sigma_0_inv, ##beta_prior variance inverse
                         block_dims, grid_dims, 
                         niter=2000, burnin=500)
{  
  beta=matrix(numeric(p*(burnin+niter+1)), (burnin+niter+1), p)
  for (i in 1: (burnin+niter))
  {
    ##rng_c related to the current iteration time
    rng_c=as.integer(i)
    
    ##calculate mean for z
    mu=as.vector(X%*%beta[i, ])
    
    ##sample z using CUDA kernel
    z= .cuda(k_rtruncnorm, "x"=val, N, mu, sigma, lo, hi, 
             rng_a, rng_b, rng_c, maxtries, m, k_rtruncnorm,
             gridDim=grid_dims, blockDim=block_dims, outputs="x")
    
    ##Calculate beta posterior mean and variance  
    po_beta_mu=solve(Sigma.0.inv+t(X)%*%X)%*%
      (Sigma.0.inv%*%beta.0+t(X)%*%z)
    po_beta_sigma=solve(Sigma.0.inv+t(X)%*%X)
    
    ##Sample beta from multivariate normal
    beta[i+1,]=mvrnorm(n = 1, mu=po_beta_mu, Sigma=po_beta_sigma)
    
    ##Print results for every 250 iterations
    if(i%%250==0)
    {
      cat("Finished", i, "Iterations...\n")
    }
  }
  cat("Work Completed...\n")
  return(beta[(burnin+2):(burnin+niter+1),])
}

##d
names=c("data_01.txt", "data_02.txt", "data_03.txt", 
        "data_04.txt", "data_05.txt")
#Compute CPU time for data_01 to data_05
compute_cpu_time=function(names)
{
  library(truncnorm)
  library(MASS)
  cat("Begin Computing: ", names, " \n")
  
  ##Set up Parameters:
  data=read.table(names, header=TRUE, sep=" ")
  y=data$y
  X=as.matrix(data[,-1])
  p=ncol(X)
  N=length(y)
  sigma=rep(1, N)
  lo=numeric(N)
  lo[which(y==0)]=-Inf
  hi=numeric(N)
  hi[which(y==1)]=Inf
  beta.0=rep(0,p)
  Sigma.0.inv=matrix(rep(0, p*p), p, p)
  
  ##Record computing time:
  time=system.time({betas=probit_mcmc_cpu(y, X, p, N, sigma, lo, hi, 
                                          beta.0, Sigma.0.inv)})
  cat("Time for ", names, " is ", time, "...\n")
  return(time)
}
cpu=sapply(names, function(x) compute_cpu_time(x))

#Compute GPU time for data_01 to data_05
compute_gpu_time=function(names)
{
  cat("Begin Computing: ", names, " \n")
  
  ##Set up Parameters:
  data=read.table(names, header=TRUE, sep=" ")  
  y=data$y
  X=as.matrix(data[,-1])
  p=ncol(X)
  N=length(y)
  sigma=rep(1, N)
  val=rep(0, N)
  lo=numeric(N)
  lo[which(y==0)]=-Inf
  hi=numeric(N)
  hi[which(y==1)]=Inf
  beta.0=rep(0,p)
  Sigma.0.inv=matrix(rep(0, p*p), p, p)
  rng_a=37L
  rng_b=46L
  maxtries=2000L
  
  ##Calculate Block, Grid dimensions
  bg = compute_grid(N)
  grid_dims = bg$grid_dims
  block_dims = bg$block_dims
  
  m = loadModule("truncnorm.ptx") 
  k_rtruncnorm = m$rtruncnorm_kernel
  
  ##Calculate running time
  time=system.time({betas=probit_mcmc_gpu(y, X, p, N, sigma, val, 
                                          lo, hi, rng_a, rng_b, maxtries, m, k_rtruncnorm,
                                          beta.0, Sigma.0.inv, block_dims, grid_dims)})
  cat("Time for ", names, " is ", time, "...\n")
  return(time)
}
gpu=sapply(names, function(x) compute_gpu_time(x))
