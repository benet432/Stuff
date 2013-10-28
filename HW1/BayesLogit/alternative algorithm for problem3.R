### Alternative Algorithm
bayes1= function(y, X, beta.star, beta.0, Sigma.0.inv)  
{
  p= ncol(X)
  pi= -0.5*matrix(beta.star-beta.0, 1, p)%*%
    solve(Sigma.0.inv)%*%as.matrix(beta.star-beta.0,p,1)+
    t(as.matrix(y))%*%data.matrix(X)%*%matrix(beta.star,p,1)-
    sum(log(1+exp(data.matrix(X)%*%matrix(beta.star,p,1))))
  return(pi)
}
bayes.logreg = function(y, X, beta.0, Sigma.0.inv, niter=10000, 
                        burnin=1000, retune=200, verbose=TRUE)
{
  p= ncol(X)
  v_2=diag(1, p)
  Sigma_new=v_2%*%Sigma.beta
  beta=matrix(numeric(p*(burnin+niter+1)), (burnin+niter+1), p)
  beta[1,]= mean_start
  U=matrix(runif(((burnin+niter)*p), 0, 1), (burnin+niter), p)
  n=rep(0, p)
  while (verbose)
  {
    for (i in 2: (retune+1))
    {      
      beta[i, ]=beta[i-1,]      
      for (j in 1:p)
      {
        beta.temp=beta[i, ]
        beta[i, j]=rnorm(1, mean=beta[i, j], sd=sqrt(Sigma_new[j, j]))
        alpha=bayes1(y, X, beta[i, ], beta.0, Sigma.0.inv)-
          bayes1(y, X, beta.temp, beta.0, Sigma.0.inv)
        if (log(U[i-1, j])<alpha)
        {
          beta[i, j]=beta[i, j] 
          n[j]=n[j]+1
        }
        else
        { beta[i, j]=beta[i-1, j] } 
      }
    }
    accept_rate=n/retune  
    n=rep(0, p)  
    
    if (all(accept_rate>=0.35) & all(accept_rate<=0.55))
    { 
      verbose = FALSE
      cat("Accept Rate for each beta is: ",
          accept_rate," We take v^2 as: ", diag(v_2), "\n")
      for (i in (2+retune): (burnin+niter+1))
      {         
        beta[i, ]=beta[i-1,]      
        for (j in 1:p)
        {
          beta.temp=beta[i, ]
          beta[i, j]=rnorm(1, mean=beta[i, j], sd=sqrt(Sigma_new[j, j]))
          alpha=bayes1(y, X, beta[i, ], beta.0, Sigma.0.inv)-
            bayes1(y, X, beta.temp, beta.0, Sigma.0.inv)
          if (log(U[i-1, j])<alpha)
          {
            beta[i, j]=beta[i, j] 
          }
          else
          { 
            beta[i, j]=beta[i-1, j] 
          } 
        }
      }              
    }
    
    else 
    {
      for (j in 1:p)
      {
        if (accept_rate[j]<0.35)
        {
          v_2[j, j]=v_2[j, j]*exp(-0.2)
        }
        else if(accept_rate[j]>0.55)
        {
          v_2[j, j]=v_2[j, j]*exp(0.2)
        }
        
      }
      cat("Accept Rate is: ",accept_rate," 
              we adjust v^2 to: ", diag(v_2), "\n")
      Sigma_new=v_2%*%Sigma.beta
    }
    
  }
  beta=beta[((burnin+2):(burnin+niter+1)), ]
  return(beta)
}

##data set-up
data=read.table(file="breast_cancer.txt", header=TRUE)
index_M=which(data[, 11]=="M")
diagnosis=numeric(nrow(data))
diagnosis[index_M]=1
data=data[, -11]
data=data.frame(data, diagnosis=diagnosis)
beta.0=matrix(rep(0, 11), 11, 1)
Sigma.0.inv=diag(1000, 11)
mean_start=glm(diagnosis~ . , data , family = "binomial")$coefficients
data=data.frame(intercept=1, data)
y=data$diagnosis
X=data[,-12]
Sigma.beta=diag(1, 11)

betas=bayes.logreg(y, X, beta.0, Sigma.0.inv)
effectiveSize(betas)