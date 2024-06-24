# testing whether new VE cpp work
rm(list=ls())
library(glmnet)
library(parallel)
numCores=detectCores()
library(Rcpp)
library(RcppEigen)
options(rcpp.warnNoExports = FALSE)
sourceCpp("estimate_VE_Eigencpp.cpp")

m=5
delta=4
p=400
n=p*delta
kappa2=0
M=m*n
kappa1=2

beta_true=rnorm(p,sd=kappa1)
X=matrix(rnorm(p*n,mean=0,sd=sqrt(1/p)),nc=p)
Xstar=matrix(rnorm(p*M,mean=0,sd=sqrt(1/p)),nc=p)

y1=rbinom(n,1, 1/(1+exp(-X%*%beta_true)))

y2=rbinom(M,1, 1/(1+exp(-Xstar%*%rep(0,p))))
y=c(y1,y2)
tau_0=0.2
wt=c(rep(1,n),rep(tau_0/m,M)) 
fit_cat<-glmnet(rbind(X,Xstar),y,weights = wt,family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)

estimate_V_ELOOCV_cpp(X,y1,Xstar,y2,tau_0,as.numeric(fit_cat$beta))

for(tau_0 in c(0.1,0.4,0.8,1.3,1.7)){
  wt=c(rep(1,n),rep(tau_0/m,M)) 
  fit_cat<-glmnet(rbind(X,Xstar),y,weights = wt,family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
  
  print(estimate_V_ELOOCV_cpp(X,y1,Xstar,y2,tau_0,as.numeric(fit_cat$beta)))
  
}



