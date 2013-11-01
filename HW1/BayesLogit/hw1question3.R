##Problem 3
## part a
##Problem 3
## part a
library(MASS)

bayes1= function(y, X, beta.star, beta.0, Sigma.0.inv) 
{
  p= ncol(X)
  pi= -0.5*matrix(beta.star-beta.0, 1, p)%*%
    solve(Sigma.0.inv)%*%as.matrix(beta.star-beta.0,p,1)+
    t(as.matrix(y))%*%data.matrix(X)%*%matrix(beta.star,p,1)-
    sum(log(1+exp(data.matrix(X)%*%matrix(beta.star,p,1))))
  return(pi)
}

bayes.logreg = function(y, X, beta.0, Sigma.0.inv, niter=50000, 
                        burnin=1000, retune=200, verbose=TRUE)
{
  p= ncol(X)
  v_2=0.01
  Sigma_new=v_2*Sigma.beta
  beta=matrix(numeric(p*(burnin+niter+1)), (burnin+niter+1), p)
  beta[1,]= mean_start
  U=runif((burnin+niter), 0, 1)
  n=0
  while (verbose)
  {
    for (i in 2: (retune+1))
    {
      
      beta[i, ]=mvrnorm(1, mu=beta[i-1, ], Sigma=Sigma_new)
      alpha=bayes1(y, X, beta[i, ], beta.0, Sigma.0.inv)-
        bayes1(y, X, beta[i-1, ], beta.0, Sigma.0.inv)
      if (log(U[i-1])<alpha)
      { beta[i, ]=beta[i, ] 
        n=n+1
      }
      else
      { beta[i, ]=beta[i-1, ] }     
    }  
    
    accept_rate=n/retune
    n=0
    
    if (accept_rate>=0.35 & accept_rate<=0.55)
    { 
      verbose = FALSE
      cat("Accept Rate is: ",accept_rate,
          " We take v^2 as: ", Sigma_new[1, 1], "\n")
      for (i in (2+retune): (burnin+niter+1))
      {
        beta[i, ]=mvrnorm(1, mu=beta[i-1, ], Sigma=Sigma_new)
        alpha=bayes1(y, X, beta[i, ], beta.0, Sigma.0.inv)-
          bayes1(y, X, beta[i-1, ], beta.0, Sigma.0.inv)
        if (log(U[i-1])<alpha)
        { beta[i, ]=beta[i, ] }
        else
        { beta[i, ]=beta[i-1, ] }       
      }  
    }
    else if (accept_rate<0.35)
    {
      verbose = TRUE
      Sigma_new=Sigma_new*exp(-0.5)
      cat("Accept Rate is: ",accept_rate,
          " we adjust v^2 to: ", Sigma_new[1, 1], "\n")
    }
    else
    {
      verbose = TRUE
      Sigma_new=Sigma_new*exp(0.5)
      cat("Accept Rate is: ",accept_rate,
          " we adjust v^2 to: ", Sigma_new[1, 1], "\n")
    }  
  }
  beta=beta[((burnin+2):(burnin+niter+1)), ]
  
  return(beta)
}

##data setup
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
Sigma.beta=summary(glm(diagnosis~ . , data , 
                       family = "binomial"))$cov.unscaled


betas=bayes.logreg(y, X, beta.0, Sigma.0.inv)
par(mfrow = c(2, 3), pty="m")
for (i in 1:11)
{
  plot(betas[,i], type="l", ylab=paste("beta_", i-1, sep=""), 
       main=paste("Traceplot for beta_", i-1, sep=""), col="blue")
}
library(coda)
effectiveSize(betas)

## part b
auto_cor=numeric(11)
for (i in 1:11)
{
  auto_cor[i]=acf(betas[,i])[1]
}
unlist(auto_cor)

## part c
cor(data)

quantiles=matrix(numeric(22), 11, 2)
for (i in 1:11)
{
  quantiles[i, ]=quantile(betas[, i], probs=c(0.025, 0.975))
}

library(BMA)
data=data[,-1]
model_selection = bic.glm(diagnosis ~ ., glm.family="binomial", data=data)
summary(model_selection)

## part d & e
p=data.matrix(X)%*%t(betas)
p=exp(p)/(1+exp(p))
p=as.vector(p)
y_simu=sapply(p, function(x) rbinom(1, 1, prob=x))
y_simu=matrix(y_simu, 569, 50000)
simu_mean=colMeans(y_simu)
simu_median=apply(y_simu, 2,  median)
ratio=sum(y==1)/sum(y==0)
simu_ratio= apply(y_simu, 2, function(x) sum(x==1)/sum(x==0))

par(mfrow = c(2, 2), pty="m")
hist(simu_mean, col="light blue", 
     main="Histogram of Simulated vs Real Data Mean", 
     xlab="Simulated Mean of Y")
abline(h=0)
abline(v=mean(y), col="red", lwd=2)

hist(simu_median, xlim=c(0,1), col="light blue", 
     main="Histogram of Simulated vs Real Data Median", 
     xlab="Simulated Median of Y")
abline(h=0)
abline(v=median(y), col="red", lwd=2)

hist(simu_ratio,  col="light blue", 
     main="Histogram of Simulated vs Real Data Ratio(1/0)", 
     xlab="Simulated Ration(1/0) of Y")
abline(h=0)
abline(v=ratio, col="red", lwd=2)