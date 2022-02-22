#
#  Multiple Linear Bayesian Regression  (iid errors) w/ DIC calculations.
#

norm.reg.mcmc <- function(y,X,beta.mean.0,beta.variance.0,n.mcmc,sigma.tune,mu.logsigma,sigma.logsigma,logsigma,n.burn,X.pred,no.print=FALSE){
  #
  # Subroutines
  #
  
  library(mvtnorm)
  
  #
  # Hyperpriors
  #
  n.pred=dim(X.pred)[1]
  n=length(y) #number of observations
  p=length(beta.mean.0) # number of predictors
  
  Sig.beta=beta.variance.0*diag(p) #covariance matrix (diagonal)
  
  beta.save=matrix(0,p,n.mcmc) #initilize beta storage
  sigmasq.save=rep(0,n.mcmc) #initilize sigmasq storage
  Dbar.save=rep(0,n.mcmc) #initialize Dbar storage
  y.pred.mean=rep(0,n) #initilize y storage
  y.pred.save=matrix(0,n.pred,n.mcmc)
  y.fit.save=matrix(0,n.pred,n.mcmc)
  
  #
  # Starting values
  #
  
  sigma = exp(logsigma) # initial value in mcmc chain
  beta=solve(t(X)%*%X)%*%t(X)%*%y # initial value (least square value)
  
  #
  # MCMC loop
  #
  
  for(k in 1:n.mcmc){
    
    # Sample log(sigma) using random walk Metropolis-Hasting 
    # then transforming samples to sigmasq
    
    logsigma.proposed = rnorm(1,logsigma,sigma.tune)
    sigma.proposed = exp(logsigma.proposed)
    #mh1 =  sigma.proposed^(-n)*exp((-1/2)*t(y-X%*%beta)%*%(y-X%*%beta)/sigma.proposed^2)*exp(-(1/2)*(1/sigma.logsigma^2)*(logsigma.proposed-mu.logsigma)^2)
    #mh2 = sigma^(-n)*exp((-1/2)*t(y-X%*%beta)%*%(y-X%*%beta)/sigma^2)*exp(-(1/2)*(1/sigma.logsigma^2)*(logsigma-mu.logsigma)^2)
    #alpha = min(1,mh1/mh2)
    #cat(alpha)
    
    mh.1 = sum(dnorm(y,X%*%beta,sigma.proposed,log=TRUE))+dnorm(logsigma.proposed,mu.logsigma,sigma.logsigma,log=TRUE)
    mh.2 = sum(dnorm(y,X%*%beta,sigma,log=TRUE))+dnorm(logsigma,mu.logsigma,sigma.logsigma,log=TRUE)
    mh = exp(mh.1-mh.2)*(sigma.proposed/sigma)^(-n)
    alpha = min(1,mh)
    if(alpha>runif(1)){
      logsigma = logsigma.proposed
    }
    sigmasq = exp(logsigma)^2
    sigma = sqrt(sigmasq)
    logsigma = log(sigma)
    
    # Sample betas from multivariate normal
    tmp.variance = solve(t(X)%*%X/sigmasq + solve(Sig.beta))
    tmp.mean = tmp.variance%*%(t(X)%*%y/sigmasq + solve(Sig.beta)%*%beta.mean.0)
    beta = as.vector(rmvnorm(1,tmp.mean,tmp.variance,method='chol'))

    # DIC calculation
    Dbar.save[k] = -2*sum(dnorm(y,X%*%beta,sqrt(sigmasq),log=TRUE))
    
    if(k>n.burn){
      y.pred = rnorm(n,X%*%beta,sqrt(sigmasq))
      y.pred.mean=y.pred.mean+y.pred/(n.mcmc-n.burn)
    }
    
    
    ###
    ### Posterior Predictive Calculations 
    ###
    
    y.fit.save[,k]=X.pred%*%beta
    y.pred.save[,k]=rnorm(n.pred,y.fit.save[,k],sigma)
    
    # Save samples
    
    beta.save[,k] = beta
    sigmasq.save[k] = sigmasq
    
    
  }
  
  #
  #  Calculate DIC and Print to Screen
  #
  
  if(dim(X)[2]==1){
    postbetamn=mean(beta.save[,-(1:n.burn)])
  }
  if(dim(X)[2]>1){
    postbetamn=apply(beta.save[,-(1:n.burn)],1,mean)
  }
  posts2mn=mean(sigmasq.save[-(1:n.burn)])
  Dhat=-2*(sum(dnorm(y,X%*%postbetamn,sqrt(posts2mn),log=TRUE)))
  Dbar=mean(Dbar.save[-(1:n.burn)])
  pD=Dbar-Dhat
  DIC=Dhat+2*pD
  
  #cat("Posterior Mean for Beta:","\n")
  #print(postbetamn)
  #cat("Posterior Mean for s2:","\n")
  #print(posts2mn)
  cat("Dhat:",Dhat,"Dbar:",Dbar,"pD:",pD,"DIC:",DIC,"\n")
  
  #
  # Write output
  #
  
  list(beta.save=beta.save,sigmasq.save=sigmasq.save,y=y,X=X,n.mcmc=n.mcmc,n=n,p=p,mu.logsigma=mu.logsigma,sigma.logsigma=sigma.logsigma,Dhat=Dhat,Dbar=Dbar,pD=pD,DIC=DIC,y.pred.mn=y.pred.mean,y.pred.save=y.pred.save,y.fit.save=y.fit.save)
}
  