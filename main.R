#
# Read data
#
setwd("~/Documents/Bayesian Statistical Methods/BayesianRegression")
car.df = read.csv("ElectricCarData_Clean.csv")
car.df

x1 = car.df[["Range_Km"]]
x2 = car.df[["Seats"]]
y = car.df[["PriceEuro"]]
n=length(y)

graphics.off()
plot(x1,y)

#
# Create design matrices
#

X.1 = model.matrix(y~x1+x2)
X.2 = model.matrix(y~x1)
X.3 = model.matrix(y~x2)

#
# Predict values
#

n.pred=100

X.1.pred=matrix(1,n.pred,3)
X.1.pred[,2]=seq(0,1000,,n.pred)

X.2.pred=matrix(1,n.pred,2)
X.2.pred[,2]=seq(0,1000,,n.pred)

X.3.pred=matrix(1,n.pred,2)
X.3.pred[,2]=seq(0,10,,n.pred)

#
# Fit model
#

sigma.tune = 0.05 # tunning parameter for random walk 
n.mcmc = 10000 # number of mcmc samples
n.burn = round(.1*n.mcmc) # burning period
logsigma = 10 # Initialize log(sigma)

#
# Hyperparameters
#

beta.mean.1.0 = c(0,0,0) # prior beta means
beta.mean.2.0 = c(0,0) # prior beta means
beta.mean.3.0 = c(0,0) # prior beta means
beta.sigma.0 = 10000000 # prior info (high variance)

#sigma.tune = 0.1 # tunning parameter for random walk 
#mu.logsigma = -1 # prior info
#sigma.logsigma = 1 # prior info

y.sigma.prior = 100000 # prior estimate of the variance of car prices given same characteristics
y.sigma.sigma.prior = 100000 # prior estimate of the variance of the variance given same characteristics
graphics.off()
data = rnorm(10000,y.sigma.prior,y.sigma.sigma.prior)
hist(data,breaks=50,prob=TRUE)

mu.logsigma = log(y.variace.prior) # prior estimate of the 
sigma.logsigma = log(y.variace.variance.prior) # prior info

mu.logsigma = log(y.sigma.prior) # prior estimate of the 
sigma.logsigma = 0.5 # prior info

graphics.off()
data = exp(rnorm(10000,mu.logsigma,sigma.logsigma))^2
hist(data,breaks=50,prob=TRUE)

source("norm.reg.mcmc.R")
mcmc.1.out = norm.reg.mcmc(y,X.1,beta.mean.1.0,beta.sigma.0,n.mcmc,sigma.tune,mu.logsigma,sigma.logsigma,logsigma,n.burn,X.1.pred,TRUE)
mcmc.2.out = norm.reg.mcmc(y,X.2,beta.mean.2.0,beta.sigma.0,n.mcmc,sigma.tune,mu.logsigma,sigma.logsigma,logsigma,n.burn,X.2.pred,TRUE)
mcmc.3.out = norm.reg.mcmc(y,X.3,beta.mean.3.0,beta.sigma.0,n.mcmc,sigma.tune,mu.logsigma,sigma.logsigma,logsigma,n.burn,X.3.pred,TRUE)

#
# Trace plots
#
layout(matrix(1:1,1,1))
plot(mcmc.1.out$beta.save[1,],type="l",lty=1,ylab=bquote(beta[0]))
plot(mcmc.1.out$beta.save[2,],type="l",lty=1,ylab=bquote(beta[1]))
plot(mcmc.1.out$beta.save[3,],type="l",lty=1,ylab=bquote(beta[2]))
plot(mcmc.1.out$sigmasq.save,type="l",ylab=bquote(sigma^2))

layout(matrix(1:1,1,1))
plot(mcmc.2.out$beta.save[1,],type="l",lty=1,ylab=bquote(beta[0]))
plot(mcmc.2.out$beta.save[2,],type="l",lty=1,ylab=bquote(beta[1]))
plot(mcmc.2.out$sigmasq.save,type="l",ylab=bquote(sigma^2))

layout(matrix(1:1,1,1))
plot(mcmc.3.out$beta.save[1,],type="l",lty=1,ylab=bquote(beta[0]))
plot(mcmc.3.out$beta.save[2,],type="l",lty=1,ylab=bquote(beta[1]))
plot(mcmc.3.out$sigmasq.save,type="l",ylab=bquote(sigma^2))

#
# Make Posterior Histograms 
#

# Km and Seats predictors
layout(matrix(1:4,2,2))
hist(mcmc.1.out$beta.save[1,],prob=TRUE,breaks=40,xlab=bquote(beta[0]),main="")
curve(dnorm(x,beta.mean.0[1],beta.sigma.0),lwd=2,add=TRUE,col=rgb(1,0,0,.5))
hist(mcmc.1.out$beta.save[2,],prob=TRUE,breaks=40,xlab=bquote(beta[1]),main="")
curve(dnorm(x,beta.mean.0[2],beta.sigma.0),lwd=2,add=TRUE,col=rgb(1,0,0,.5))
hist(mcmc.1.out$beta.save[3,],prob=TRUE,breaks=40,xlab=bquote(beta[2]),main="")
curve(dnorm(x,beta.mean.0[3],beta.sigma.0),lwd=2,add=TRUE,col=rgb(1,0,0,.5))
hist(mcmc.1.out$sigmasq.save,prob=TRUE,breaks=40,xlab=bquote(sigma^2),main="")

# Km predictor

layout(matrix(1:4,2,2))
hist(mcmc.2.out$beta.save[1,],prob=TRUE,breaks=40,xlab=bquote(beta[0]),main="")
curve(dnorm(x,beta.mean.0[1],beta.sigma.0),lwd=2,add=TRUE,col=rgb(1,0,0,.5))
hist(mcmc.2.out$beta.save[2,],prob=TRUE,breaks=40,xlab=bquote(beta[1]),main="")
curve(dnorm(x,beta.mean.0[2],beta.sigma.0),lwd=2,add=TRUE,col=rgb(1,0,0,.5))
hist(mcmc.2.out$sigmasq.save,prob=TRUE,breaks=40,xlab=bquote(sigma^2),main="")
curve(dnorm(x,mu.logsigma,sigma.logsigma),lwd=2,add=TRUE,col=rgb(1,0,0,.5))

pred.mn=apply(mcmc.2.out$y.pred.save,1,mean)
pred.CI=t(apply(mcmc.2.out$y.pred.save,1,quantile,c(0.025,0.975)))
fit.mn=apply(mcmc.2.out$y.fit.save,1,mean)
fit.CI=t(apply(mcmc.2.out$y.fit.save,1,quantile,c(0.025,0.975)))

layout(matrix(c(1,1,1,1,2,3),2,3))
matplot(X.2.pred[,2],pred.CI,type="l",col=2,lwd=2,lty=1,xlab="Range(Km)",ylab="price")
matplot(X.2.pred[,2],fit.CI,type="l",col=3,lwd=2,lty=1,add=TRUE)
points(x1,y,add=TRUE)

# Seats predictor
layout(matrix(1:4,2,2))
hist(mcmc.3.out$beta.save[1,],prob=TRUE,breaks=40,xlab=bquote(beta[0]),main="")
curve(dnorm(x,beta.mean.0[1],beta.variance.0),lwd=2,add=TRUE,col=rgb(1,0,0,.5))
hist(mcmc.3.out$beta.save[2,],prob=TRUE,breaks=40,xlab=bquote(beta[1]),main="")
curve(dnorm(x,beta.mean.0[2],beta.variance.0),lwd=2,add=TRUE,col=rgb(1,0,0,.5))
hist(mcmc.3.out$sigmasq.save,prob=TRUE,breaks=40,xlab=bquote(sigma^2),main="")


pred.mn=apply(mcmc.3.out$y.pred.save,1,mean)
pred.CI=t(apply(mcmc.3.out$y.pred.save,1,quantile,c(0.025,0.975)))
fit.mn=apply(mcmc.3.out$y.fit.save,1,mean)
fit.CI=t(apply(mcmc.3.out$y.fit.save,1,quantile,c(0.025,0.975)))

layout(matrix(c(1,1,1,1,2,3),2,3))
matplot(X.3.pred[,2],pred.CI,type="l",col=2,lwd=2,lty=1,xlab="Seats",ylab="price")
matplot(X.3.pred[,2],fit.CI,type="l",col=3,lwd=2,lty=1,add=TRUE)
points(x2,y,add=TRUE)

