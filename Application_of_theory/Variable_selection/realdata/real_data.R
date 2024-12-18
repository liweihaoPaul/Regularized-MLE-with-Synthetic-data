
gene_name<-read.table("top600_gene_name.txt",header=F)
gene_name<-as.vector(gene_name$V1)
library(MASS)
library(glmnet)
library(parallel)
num_cores <- detectCores()
library(Rcpp)
library(RcppEigen)
sourceCpp("estimate_eta_square_Eigencpp.cpp")

### source code
source("ABHqMAP.R")
source("ABYqMAP.R")
source("DSMAP.R")
source("MDSMAP.R")
source("utilsMAP.R")


f <- function(x){
  exp(x)/(1 + exp(x))
}

p=ncol(final_data)
n=nrow(final_data)
M=20*p
delta=n/p
set.seed(1)

y=c(rep(0,400),rep(1,2000))
X=final_data/sqrt(n)
X=as.matrix(X)
# generate synthetic data
sigma_hat=crossprod(X)/n
X.syn = mvrnorm(M, mu = rep(0, p), Sigma = sigma_hat)
y.syn = rbinom(M, 1, rep(0.5,M))

# find sample covariance of X, EX is column mean sigma_hat=(X-EX)^T(X-EX)/n
EX=colMeans(X)
sigma_hat=cov(X)
X.syn = mvrnorm(M, mu = EX, Sigma = sigma_hat)
y.syn = rbinom(M, 1, rep(0.5,M))


q=0.1
ABHq_time = system.time(ABHq_result <- ABHqMAP(X, y,X.syn,y.syn))[3]
ABHq_result$select_index
gene_name[ABHq_result$select_index]

ABYq_time = system.time(ABYq_result <- ABYqMAP(X, y,X.syn,y.syn))[3]
ABYq_result$select_index
gene_name[ABYq_result$select_index]

MDS_time = system.time(MDS_result <- MDSMAP(X, y, num_split = 50))[3]
MDS_result$select_index
gene_name[MDS_result$select_index]


