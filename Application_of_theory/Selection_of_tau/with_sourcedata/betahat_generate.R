# compute betahat for 50 replications, later can be used to compute mse and pred deviance
# four settings: kappa1={1,2},delta={2,4}, kesee=0.9, delta_s=10,kappa_2=1

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

betahat_LOOCV_informative_syn_p400<-function(kappa1,delta,kesee=0.9,
                                                 kappa2=1,delta_s=10,beta_true_gen="t3",N_rep=50,p_seq=c(100 , 400 ,1600)){
  repeat_function_for_betahat<-function(rep_i){
    ####### generate target data first
    set.seed(rep_i+100) # for observed data beta
    beta_true_large=generate_beta(p_seq[3],kappa1,beta_true_gen)
    set.seed(rep_i+200) # for observed data X
    X_large=matrix(rnorm(p_seq[3]^2*delta , mean = 0, sd = 1), nc = p_seq[3])
    # only for p400
    pindex=2
    set.seed(rep_i+1000*pindex)
    p=p_seq[pindex]
    n=p*delta
    
    beta_true=beta_true_large[1:p]/sqrt(p)
    
    X = X_large[1:n,1:p]
    Y = rbinom(n, 1, 1 / (1 + exp(-X %*% beta_true)))
    target_X=X
    target_Y=Y
    ####### generate source data
    set.seed(rep_i+300) # for source data beta and X
    lambda=kappa2*(sqrt(1-kesee^2))
    betas_large=kesee*kappa2/kappa1*beta_true_large+lambda*generate_beta(p_seq[3],kappa1=1,beta_true_gen)
    X_large=matrix(rnorm(p_seq[3]^2*delta_s , mean = 0, sd = 1), nc = p_seq[3])
    # only for p400
    pindex=2
    set.seed(rep_i+1000*pindex)
    p=p_seq[pindex]
    n=p*delta_s
    
    betas=betas_large[1:p]/sqrt(p)
    
    X = X_large[1:n,1:p]
    Y = rbinom(n, 1, 1 / (1 + exp(-X %*% betas)))
    source_X=X
    source_Y=Y
    best_tau_0<-loocv_best_tau(target_X,target_Y,source_X,source_Y,tau_0_seq)
    fit_cat_besttau_0_loocvsource=glmnet(rbind(target_X,source_X),c(target_Y,source_Y),
                                   weights =c(rep(1,length(target_Y)),rep(best_tau_0*length(target_Y)/length(source_Y),length(source_Y))) ,
                                   family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
    return(c(best_tau_0,as.numeric(fit_cat_besttau_0_loocvsource$beta))) # first entries is tau_0_best
  }
  parallel_betahatsourceloocv_list=mclapply(1:N_rep,repeat_function_for_betahat,mc.cores=numCores)
  save(parallel_betahatsourceloocv_list,file=paste0("betahat_LOOCV_informative_syn_p400_kappa1_",kappa1,"_delta_",delta,"_kesee_",kesee,"_kappa2_",kappa2,"_delta_s_",delta_s,".RData"))
  
  
}
  

for(kappa1 in c(1,2)){
  for(delta in c(2,4)){
    betahat_LOOCV_informative_syn_p400(kappa1,delta)
  }
}



# best tau from kappa1 hat ------------------------------------------------
betahat_estkappa1kappa2kesee_mse_curve_informative_syn_p400<-function(kappa1,delta,kesee=0.9,
                                                                      kappa2=1,delta_s=10,beta_true_gen="t3",N_rep=50,p_seq=c(100 , 400 ,1600)){
  load(paste0("result_list_solution_kappa12keseehat_p400_kappa1_",
              kappa1,"_delta_",delta,"_tau_0_",
              1/4,"_kesee_",kesee,"_kappa2_",
              kappa2,"_delta_s_",delta_s,"_beta_true_gen_",
              beta_true_gen,".RData"))
  solution_list_small_tau=result_list_solution_kappa12keseehat_p400
  load(paste0("appendresult_list_solution_kappa12keseehat_p400_kappa1_",
              kappa1,"_delta_",delta,"_tau_0_",
              1/4,"_kesee_",kesee,"_kappa2_",
              kappa2,"_delta_s_",delta_s,"_beta_true_gen_",
              beta_true_gen,".RData"))
  solution_list_large_tau=result_list_solution_kappa12keseehat_p400
  load(paste0("result_matrix_kappa2_kappa1_is_",kappa1,"_delta_",delta,"_tau_0_",
              1/4,"_kesee_",kesee,"_kappa2_",
              kappa2,"_delta_s_",delta_s,"_beta_true_gen_",
              beta_true_gen,".RData"))
  load(paste0("result_matrix_kappa1_",kappa1,
              "_delta_",delta,"_tau_0_",
              1/4,"_beta_true_gen_",beta_true_gen,".RData")) 
  kappa1_hat_seq=est_kappa1_matrix[,2]
  kappa2_hat_seq=est_kappa2_matrix[,2]
  
  repeat_function_for_betahat<-function(rep_i){
    ####### generate target data first
    set.seed(rep_i+100) # for observed data beta
    beta_true_large=generate_beta(p_seq[3],kappa1,beta_true_gen)
    set.seed(rep_i+200) # for observed data X
    X_large=matrix(rnorm(p_seq[3]^2*delta , mean = 0, sd = 1), nc = p_seq[3])
    # only for p400
    pindex=2
    set.seed(rep_i+1000*pindex)
    p=p_seq[pindex]
    n=p*delta
    
    beta_true=beta_true_large[1:p]/sqrt(p)
    
    X = X_large[1:n,1:p]
    Y = rbinom(n, 1, 1 / (1 + exp(-X %*% beta_true)))
    target_X=X
    target_Y=Y
    ####### generate source data
    set.seed(rep_i+300) # for source data beta and X
    lambda=kappa2*(sqrt(1-kesee^2))
    betas_large=kesee*kappa2/kappa1*beta_true_large+lambda*generate_beta(p_seq[3],kappa1=1,beta_true_gen)
    X_large=matrix(rnorm(p_seq[3]^2*delta_s , mean = 0, sd = 1), nc = p_seq[3])
    # only for p400
    pindex=2
    set.seed(rep_i+1000*pindex)
    p=p_seq[pindex]
    n=p*delta_s
    
    betas=betas_large[1:p]/sqrt(p)
    
    X = X_large[1:n,1:p]
    Y = rbinom(n, 1, 1 / (1 + exp(-X %*% betas)))
    source_X=X
    source_Y=Y
    
    # compute best_tau_0
    kappa1_hat=kappa1_hat_seq[rep_i]
    kappa2_hat=kappa2_hat_seq[rep_i]
    TheorySolution1=rbind(solution_list_small_tau[[rep_i]],solution_list_large_tau[[rep_i]])
    mse_cirve_estkappa1kappa2kesee_seq=TheorySolution1[,1]^2+kappa1_hat^2*(TheorySolution1[,3]-1)^2+kappa2_hat^2*TheorySolution1[,4]^2
    best_tau_0=tau_0_seq[which.min(mse_cirve_estkappa1kappa2kesee_seq)]
    
    fit_cat_besttau_0_estkesee=glmnet(rbind(target_X,source_X),c(target_Y,source_Y),
                                   weights =c(rep(1,length(target_Y)),rep(best_tau_0*length(target_Y)/length(source_Y),length(source_Y))) ,
                                   family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
    return(c(best_tau_0,as.numeric(fit_cat_besttau_0_estkesee$beta))) # first entries is tau_0_best
  }
  parallel_betahatsourceestkesee_mse_list=mclapply(1:N_rep,repeat_function_for_betahat,mc.cores=numCores)
  save(parallel_betahatsourceestkesee_mse_list,file=paste0("betahat_estkappa1kappa2kesee_mse_curve_informative_syn_p400_kappa1_",kappa1,"_delta_",delta,"_kesee_",kesee,"_kappa2_",kappa2,"_delta_s_",delta_s,".RData"))
}

for(kappa1 in c(1,2)){
  for(delta in c(2,4)){
    betahat_estkappa1kappa2kesee_mse_curve_informative_syn_p400(kappa1,delta)
  }
}


# best tau from oracle kappa1 kesee known ------------------------------------------------
betahat_oraclekesee_mse_curve_informative_syn_p400<-function(kappa1,delta,kesee=0.9,
                                                                      kappa2=1,delta_s=10,beta_true_gen="t3",N_rep=50,p_seq=c(100 , 400 ,1600)){
  load(paste0("TheorySolution_infor_source_data_delta_kappa1_",
              kappa1,"_kappa2_",kappa2,"_delta_",delta,"_kesee_",kesee,"_delta_s_",10,".RData"))
  TheorySolution_small_tau=TheorySolution
  load(paste0("appendTheorySolution_infor_source_data_delta_kappa1_",
              kappa1,"_kappa2_",kappa2,"_delta_",delta,"_kesee_",kesee,"_delta_s_",10,".RData"))
  TheorySolution_large_tau=TheorySolution
  TheorySolution1=rbind(TheorySolution_small_tau,TheorySolution_large_tau)
  mse_cirve_oraclekesee_seq=TheorySolution1[,1]^2+kappa1^2*(TheorySolution1[,3]-1)^2+kappa2^2*TheorySolution1[,4]^2
  best_tau_0=tau_0_seq[which.min(mse_cirve_oraclekesee_seq)]
  
  repeat_function_for_betahat<-function(rep_i){
    ####### generate target data first
    set.seed(rep_i+100) # for observed data beta
    beta_true_large=generate_beta(p_seq[3],kappa1,beta_true_gen)
    set.seed(rep_i+200) # for observed data X
    X_large=matrix(rnorm(p_seq[3]^2*delta , mean = 0, sd = 1), nc = p_seq[3])
    # only for p400
    pindex=2
    set.seed(rep_i+1000*pindex)
    p=p_seq[pindex]
    n=p*delta
    
    beta_true=beta_true_large[1:p]/sqrt(p)
    
    X = X_large[1:n,1:p]
    Y = rbinom(n, 1, 1 / (1 + exp(-X %*% beta_true)))
    target_X=X
    target_Y=Y
    ####### generate source data
    set.seed(rep_i+300) # for source data beta and X
    lambda=kappa2*(sqrt(1-kesee^2))
    betas_large=kesee*kappa2/kappa1*beta_true_large+lambda*generate_beta(p_seq[3],kappa1=1,beta_true_gen)
    X_large=matrix(rnorm(p_seq[3]^2*delta_s , mean = 0, sd = 1), nc = p_seq[3])
    # only for p400
    pindex=2
    set.seed(rep_i+1000*pindex)
    p=p_seq[pindex]
    n=p*delta_s
    
    betas=betas_large[1:p]/sqrt(p)
    
    X = X_large[1:n,1:p]
    Y = rbinom(n, 1, 1 / (1 + exp(-X %*% betas)))
    source_X=X
    source_Y=Y
    
    fit_cat_besttau_0_oraclekesee=glmnet(rbind(target_X,source_X),c(target_Y,source_Y),
                                      weights =c(rep(1,length(target_Y)),rep(best_tau_0*length(target_Y)/length(source_Y),length(source_Y))) ,
                                      family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
    return(c(best_tau_0,as.numeric(fit_cat_besttau_0_oraclekesee$beta))) # first entries is tau_0_best
  }
  parallel_betahatsourceoraclekesee_mse_list=mclapply(1:N_rep,repeat_function_for_betahat,mc.cores=numCores)
  save(parallel_betahatsourceoraclekesee_mse_list,file=paste0("betahat_oraclekesee_mse_curve_informative_syn_p400_kappa1_",
                                                              kappa1,"_delta_",delta,"_kesee_",kesee,"_kappa2_",kappa2,"_delta_s_",delta_s,".RData"))
}

for(kappa1 in c(1,2)){
  for(delta in c(2,4)){
    betahat_oraclekesee_mse_curve_informative_syn_p400(kappa1,delta)
  }
}








