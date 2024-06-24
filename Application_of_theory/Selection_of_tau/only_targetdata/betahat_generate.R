# compute betahat for 50 replications, later can be used to compute mse and pred deviance
# four settings: kappa1={1,2},delta={2,4}
library(glmnet)
library(parallel)
numCores=min(50,detectCores())
library(Rcpp)
library(RcppEigen)
options(rcpp.warnNoExports = FALSE)
sourceCpp("estimate_VE_Eigencpp.cpp")

tau_0_seq=seq(0.1,4,0.1)
generate_beta<-function(p,kappa1,beta_true_gen="normal"){
  if(beta_true_gen=="normal"){
    return(rnorm(p,sd=kappa1))
  }
  if(beta_true_gen=="uniform"){
    return(runif(p,min=-1,max=1)*kappa1*sqrt(3))
  }
  if(beta_true_gen=="t3"){
    return(rt(p,df=3)*kappa1/sqrt(3))
  }
  if(beta_true_gen=="halfsparse"){
    ind <- sample(c(TRUE, FALSE), p, replace = TRUE)
    v=rnorm(p,sd=kappa1*sqrt(2))
    v[ind]=0
    return(v)
  }
}

# leave one out cross validation ------------------------------------------

loocv_best_tau<-function(X,Y,Xstar,Ystar,tau_0_seq){
  p=ncol(X)
  n=length(Y)
  M=length(Ystar)
  ve_error=rep(0,length(tau_0_seq))
  for(tau_0_index in 1:length(tau_0_seq)){
    tau_0=tau_0_seq[tau_0_index]
    fit_cat_tau_0=glmnet(rbind(X,Xstar),c(Y,Ystar),weights =c(rep(1,n),rep(tau_0*n/M,M)) ,
                         family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
    betahat_tau0=as.numeric(fit_cat_tau_0$beta)
    ve_error[tau_0_index]=estimate_VELOOCV_cpp(X,Y,Xstar,Ystar,tau_0,betahat_tau0)
  }
  return(tau_0_seq[which.min(ve_error)])
}

# this function will save a list with length N_rep, each element is a vector of betahat
betahat_LOOCV_non_informative_syn_p400<-function(kappa1,delta,beta_true_gen="t3",N_rep=50,M_div_p=20,p_seq=c(100 , 400 ,1600)){
  repeat_function_for_betahat<-function(rep_i){
    set.seed(rep_i+100) # for observed data beta
    beta_true_large=generate_beta(p_seq[3],kappa1,beta_true_gen)
    set.seed(rep_i+200) # for observed data X
    X_large=matrix(rnorm(p_seq[3]^2*delta , mean = 0, sd = 1), nc = p_seq[3])
    # set.seed(rep_i+300) # for source data beta and X
    # betas_large=kesee*kappa2/kappa1*beta_true_large+lambda*generate_beta(p_seq[3],kappa1=1,beta_true_gen)
    # X_large=matrix(rnorm(p_seq[3]^2*delta_s , mean = 0, sd = 1), nc = p_seq[3]) # make sure generate different data
    set.seed(rep_i+400) # for synthetic data
    Xstar_large=matrix(rnorm(p_seq[3]^2*M_div_p , mean = 0, sd = 1), nc = p_seq[3])
    # only for p400
    pindex=2
    set.seed(rep_i+1000*pindex)
    p=p_seq[pindex]
    n=p*delta
    M=M_div_p*p
    beta_true=beta_true_large[1:p]/sqrt(p)
    
    X = X_large[1:n,1:p]
    Y = rbinom(n, 1, 1 / (1 + exp(-X %*% beta_true)))
    Xstar = Xstar_large[1:M,1:p]
    beta0 = rep(0, p)
    Ystar = rbinom(M, 1, 1 / (1 + exp(-Xstar %*% beta0)))
    best_tau_0<-loocv_best_tau(X,Y,Xstar,Ystar,tau_0_seq)
    fit_cat_besttau_0_loocv=glmnet(rbind(X,Xstar),c(Y,Ystar),weights =c(rep(1,n),rep(best_tau_0*n/M,M)) ,
                                   family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
    return(c(best_tau_0,as.numeric(fit_cat_besttau_0_loocv$beta))) # first entries is tau_0_best
  }
  parallel_betahatloocv_list=mclapply(1:N_rep,repeat_function_for_betahat,mc.cores=numCores)
  save(parallel_betahatloocv_list,file=paste0("parallel_betahatloocv_list_non_informative_syn_p400_kappa1_",kappa1,"_delta_",delta,".RData"))
}

for(kappa1 in c(1,2)){
  for(delta in c(2,4)){
    betahat_LOOCV_non_informative_syn_p400(kappa1,delta)
  }
}


# best tau from kappa1 hat ------------------------------------------------


betahat_estkappa1_mse_curve_non_informative_syn_p400<-function(kappa1,delta,beta_true_gen="t3",N_rep=50,M_div_p=20,p_seq=c(100 , 400 ,1600)){
  load(paste0("result_list_solution_kappahat_p400_kappa1_",kappa1,
              "_delta_",delta,".RData"))
  solutionlist_smalltau=result_list_solution_kappahat_p400
  load(paste0("appendresult_list_solution_kappahat_p400_kappa1_",kappa1,
              "_delta_",delta,".RData"))
  solutionlist_largetau=result_list_solution_kappahat_p400
  load(paste0("result_matrix_kappa1_",kappa1,
                   "_delta_",delta,"_tau_0_",
                   1/4,"_beta_true_gen_","t3",".RData"))
  
  
  kappa1_hat_seq=est_kappa1_matrix[,2]
  repeat_function_for_betahat<-function(rep_i){
    kappa1_hat=kappa1_hat_seq[rep_i]
    solution_matrix_different_tau_0=rbind(solutionlist_smalltau[[rep_i]],solutionlist_largetau[[rep_i]])
    msehat_tau_0_seq=solution_matrix_different_tau_0[,1]^2+kappa1_hat^2*(solution_matrix_different_tau_0[,3]-1)^2
    best_tau_0=tau_0_seq[which.min(msehat_tau_0_seq)]
    
    set.seed(rep_i+100) # for observed data beta
    beta_true_large=generate_beta(p_seq[3],kappa1,beta_true_gen)
    set.seed(rep_i+200) # for observed data X
    X_large=matrix(rnorm(p_seq[3]^2*delta , mean = 0, sd = 1), nc = p_seq[3])
    # set.seed(rep_i+300) # for source data beta and X
    # betas_large=kesee*kappa2/kappa1*beta_true_large+lambda*generate_beta(p_seq[3],kappa1=1,beta_true_gen)
    # X_large=matrix(rnorm(p_seq[3]^2*delta_s , mean = 0, sd = 1), nc = p_seq[3]) # make sure generate different data
    set.seed(rep_i+400) # for synthetic data
    Xstar_large=matrix(rnorm(p_seq[3]^2*M_div_p , mean = 0, sd = 1), nc = p_seq[3])
    # only for p400
    pindex=2
    set.seed(rep_i+1000*pindex)
    p=p_seq[pindex]
    n=p*delta
    M=M_div_p*p
    beta_true=beta_true_large[1:p]/sqrt(p)
    
    X = X_large[1:n,1:p]
    Y = rbinom(n, 1, 1 / (1 + exp(-X %*% beta_true)))
    Xstar = Xstar_large[1:M,1:p]
    beta0 = rep(0, p)
    Ystar = rbinom(M, 1, 1 / (1 + exp(-Xstar %*% beta0)))
    fit_cat_besttau_0_estkappa1_mse_curve=glmnet(rbind(X,Xstar),c(Y,Ystar),weights =c(rep(1,n),rep(best_tau_0*n/M,M)) ,
                                   family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
    return(c(best_tau_0,as.numeric(fit_cat_besttau_0_estkappa1_mse_curve$beta))) # first entries is tau_0_best
  }
  parallel_betahatestkappa1_mse_list=mclapply(1:N_rep,repeat_function_for_betahat,mc.cores=numCores)
  save(parallel_betahatestkappa1_mse_list,file=paste0("parallel_betahatestkappa1_mse_list_non_informative_syn_p400_kappa1_",kappa1,"_delta_",delta,".RData"))
}

for(kappa1 in c(1,2)){
  for(delta in c(2,4)){
    betahat_estkappa1_mse_curve_non_informative_syn_p400(kappa1,delta)
  }
}


# oracle tau from mse curve -----------------------------------------------
betahat_oracle_mse_curve_non_informative_syn_p400<-function(kappa1,delta,beta_true_gen="t3",N_rep=50,M_div_p=20,p_seq=c(100 , 400 ,1600)){
  load(paste0("TheorySolution_non_infor_syn_data_delta_kappa1_",kappa1,"_delta_",delta,".RData"))
  Theorysolution_smalltau=TheorySolution
  load(paste0("appendTheorySolution_non_infor_syn_data_delta_kappa1_",kappa1,"_delta_",delta,".RData"))
  Theorysolution_largetau=TheorySolution
  TheorySolution=rbind(Theorysolution_smalltau,Theorysolution_largetau)
  mse_seq_diff_tau_0=TheorySolution[,1]^2+kappa1^2*(TheorySolution[,3]-1)^2
  best_tau_0=tau_0_seq[which.min(mse_seq_diff_tau_0)]
  repeat_function_for_betahat<-function(rep_i){
    
    set.seed(rep_i+100) # for observed data beta
    beta_true_large=generate_beta(p_seq[3],kappa1,beta_true_gen)
    set.seed(rep_i+200) # for observed data X
    X_large=matrix(rnorm(p_seq[3]^2*delta , mean = 0, sd = 1), nc = p_seq[3])
    # set.seed(rep_i+300) # for source data beta and X
    # betas_large=kesee*kappa2/kappa1*beta_true_large+lambda*generate_beta(p_seq[3],kappa1=1,beta_true_gen)
    # X_large=matrix(rnorm(p_seq[3]^2*delta_s , mean = 0, sd = 1), nc = p_seq[3]) # make sure generate different data
    set.seed(rep_i+400) # for synthetic data
    Xstar_large=matrix(rnorm(p_seq[3]^2*M_div_p , mean = 0, sd = 1), nc = p_seq[3])
    # only for p400
    pindex=2
    set.seed(rep_i+1000*pindex)
    p=p_seq[pindex]
    n=p*delta
    M=M_div_p*p
    beta_true=beta_true_large[1:p]/sqrt(p)
    
    X = X_large[1:n,1:p]
    Y = rbinom(n, 1, 1 / (1 + exp(-X %*% beta_true)))
    Xstar = Xstar_large[1:M,1:p]
    beta0 = rep(0, p)
    Ystar = rbinom(M, 1, 1 / (1 + exp(-Xstar %*% beta0)))
    fit_cat_besttau_0_oraclemse=glmnet(rbind(X,Xstar),c(Y,Ystar),weights =c(rep(1,n),rep(best_tau_0*n/M,M)) ,
                                   family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
    return(c(best_tau_0,as.numeric(fit_cat_besttau_0_oraclemse$beta))) # first entries is tau_0_best
  }
  parallel_betahateoracle_mse_list=mclapply(1:N_rep,repeat_function_for_betahat,mc.cores=numCores)
  save(parallel_betahateoracle_mse_list,file=paste0("parallel_betahateoracle_mse_list_non_informative_syn_p400_kappa1_",kappa1,"_delta_",delta,".RData"))
}

for(kappa1 in c(1,2)){
  for(delta in c(2,4)){
    betahat_oracle_mse_curve_non_informative_syn_p400(kappa1,delta)
  }
}

