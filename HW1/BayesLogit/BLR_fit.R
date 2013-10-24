##
#
# Logistic regression
# 
# Y_{i} | \beta \sim \textrm{Bin}\left(n_{i},e^{x_{i}^{T}\beta}/(1+e^{x_{i}^{T}\beta})\right)
# \beta \sim N\left(\beta_{0},\Sigma_{0}\right)
#
##
library(MASS)
library(mvtnorm)
library(coda)

########################################################################################
########################################################################################
## Handle batch job arguments:

# 1-indexed version is used now.
args <- commandArgs(TRUE)

cat(paste0("Command-line arguments:\n"))
print(args)

####
# sim_start ==> Lowest simulation number to be analyzed by this particular batch job
###

#######################
sim_start <- 1000
length.datasets <- 200
#######################

if (length(args)==0){
  sinkit <- FALSE
  sim_num <- sim_start + 1
  set.seed(1330931)
} else {
  # Sink output to file?
  sinkit <- TRUE
  # Decide on the job number, usually start at 1000:
  sim_num <- sim_start + as.numeric(args[1])
  # Set a different random seed for every job number!!!
  set.seed(762*sim_num + 1330931)
}

# Simulation datasets numbered 1001-1200

########################################################################################
########################################################################################

"bayes1"<- function(m, y, X, beta.star, beta.0, Sigma.0.inv)
  
{
  n= length(m)
  p= ncol(X)
  pi= -0.5*matrix(beta.star-beta.0, 1, p)%*%Sigma.0.inv%*%as.matrix(beta.star-beta.0,p,1)+matrix(y, 1, n)%*%X%*%matrix(beta.star,p,1)-matrix(m, 1, n)%*%log(1+exp(X%*%matrix(beta.star,p,1)))
  return(pi)
}

"bayes.logreg" <- function(m, y, X, beta.0=c(0, 0), Sigma.0.inv=diag(1, 2), niter=10000, burnin=1000, print.every=1000, retune=200, verbose=TRUE)
{
  beta_ci=matrix(numeric(2*99), 99, 2)
  v_2=0.01
  Sigma_new=v_2*diag(1, 2)
  beta=matrix(numeric(2*(burnin+niter+1)), (burnin+niter+1), 2)
  beta[1,]= beta.0
  n= length(m)
  p= ncol(X)
  U=runif((burnin+niter), 0, 1)
  n=0
  while (verbose)
  {
    for (i in 2: (retune+1))
    {
      
      beta[i, ]=mvrnorm(1, mu=beta[i-1, ], Sigma=Sigma_new)
      alpha=bayes1(m, y, X, beta[i, ], beta.0, Sigma.0.inv)-bayes1(m, y, X, beta[i-1, ], beta.0, Sigma.0.inv)
      if (log(U[i-1])<alpha)
      { beta[i, ]=beta[i, ] 
        n=n+1
      }
      else
      { beta[i, ]=beta[i-1, ] }     
    }  
    
    accept_rate=n/retune
    n=0
    
    if (accept_rate>=0.3 & accept_rate<=0.6)
    { 
      verbose = FALSE
      cat("Accept Rate is: ",accept_rate," We take v^2 as: ", Sigma_new[1, 1], "\n")
      for (i in (2+retune): (burnin+niter+1))
      {
        beta[i, ]=mvrnorm(1, mu=beta[i-1, ], Sigma=Sigma_new)
        alpha=bayes1(m, y, X, beta[i, ], beta.0, Sigma.0.inv)-bayes1(m, y, X, beta[i-1, ], beta.0, Sigma.0.inv)
        if (log(U[i-1])<alpha)
        { beta[i, ]=beta[i, ] }
        else
        { beta[i, ]=beta[i-1, ] }       
      }  
    }
    else if (accept_rate<0.3)
    {
      verbose = TRUE
      Sigma_new=Sigma_new*exp(-0.5)
      cat(cat("Accept Rate is: ",accept_rate," we adjust v^2 to: ", Sigma_new[1, 1], "\n"))
    }
    else
    {
      verbose = TRUE
      Sigma_new=Sigma_new*exp(0.5)
      cat(cat("Accept Rate is: ",accept_rate," we adjust v^2 to: ", Sigma_new[1, 1], "\n"))
    }  
  }
  beta=beta[((burnin+2):(burnin+niter+1)), ]
  for (j in 1:2)
  {
    beta_ci[, j]=quantile(beta[, j], probs=seq(0.01, 0.99, 0.01))
  }
  return(beta_ci)
}

#################################################
# Set up the specifications:

data=read.csv(file=paste('data/blr_data_', as.character(sim_num), '.csv', sep=""), header=FALSE, sep=",")
m=data$n
y=data$y
X=as.matrix(data[, 3:4])
beta.0 = matrix(c(0,0))
Sigma.0.inv = diag(rep(1.0,p))
niter = 10000
beta_cre=bayes.logreg(m, y, X, beta.0, Sigma.0.inv=diag(1, 2), niter=10000, burnin=1000, print.every=1000, retune=200, verbose=TRUE)
write.table(beta_cre, file=paste('results/blr_res_', as.character(sim_num), '.csv', sep=""), sep=",", row.names = FALSE, col.names = FALSE)

# etc... (more needed here)
#################################################

# Read data corresponding to appropriate sim_num:

# Extract X and y:

# Fit the Bayesian model:

# Extract posterior quantiles...

# Write results to a (99 x p) csv file...

# Go celebrate.

cat("done. :)\n")