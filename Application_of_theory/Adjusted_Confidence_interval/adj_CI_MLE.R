# based on code from 'glmhd' qian zhao 2022


library(cubature)
library(ECOSolveR)
library(pracma)
library(stats)
library(Matrix)
library(robustbase)





library(parallel)

numCores=detectCores()
r_files <- list.files(pattern = "\\.R$", full.names = TRUE)
r_files<-r_files[-which(r_files=="./adj_CI_MLE.R")]

for (file in r_files) {
  source(file)
}




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




coverage_adjMLE<-function(kappa1,delta,beta_true_gen="t3",N_rep=100,M_div_p=20,p_seq=c(100 , 400 ,1600)){
  result_matrix_mlecoverage=matrix(0,nc=length(p_seq)-1,nr=N_rep)
  
  repeat_function_for_eta_square<-function(rep_i){
    print(rep_i)
    cover_result_all_p=rep(0,length(p_seq)-1)
    set.seed(rep_i+100) # for observed data beta
    beta_true_large=generate_beta(p_seq[3],kappa1,beta_true_gen)
    set.seed(rep_i+200) # for observed data X
    X_large=matrix(rnorm(p_seq[3]^2*delta , mean = 0, sd = 1), nc = p_seq[3])
    
    for(pindex in 1:(length(p_seq)-1)){
      set.seed(rep_i+1000*pindex)
      p=p_seq[pindex]
      n=p*delta
      M=M_div_p*p
      beta_true=beta_true_large[1:p]
      
      X = X_large[1:n,1:p]/sqrt(p)
      X=X/sqrt(delta)
      beta_true=beta_true*sqrt(delta)
      
      Y = rbinom(n, 1, 1 / (1 + exp(-X %*% beta_true)))
      have_warning=FALSE
      tryCatch({
        fit_glm=glm(Y~X-1,family = binomial(link = "logit"),x=TRUE,y=TRUE)
      }, warning = function(w) {
        have_warning<<-TRUE
      })
      if(have_warning==FALSE){
        adjusted_fit <- adjust_glm(fit_glm, verbose = FALSE, echo = TRUE)
        est=summary(adjusted_fit)[["coefficients"]][,1]
        sdd=summary(adjusted_fit)[["coefficients"]][,1]
        upper_CI=as.numeric(est+qnorm(0.975)*sdd)
        lower_CI=as.numeric(est-qnorm(0.975)*sdd)
        cover_result_all_p[pindex]=mean((beta_true>lower_CI)&(beta_true<upper_CI))
      }else{
        cover_result_all_p[pindex]=999
      }
    }
    return(cover_result_all_p)
  }
  
  parallel_eta_list=mclapply(1:N_rep,repeat_function_for_eta_square,mc.cores=numCores)
  for(rep_i in 1:N_rep){
    results_etasquare <- parallel_eta_list[[rep_i]]
    result_matrix_adjmlecoverage[rep_i,]=results_etasquare
  }
  
  
  # save file with kappa1 delta and tau_0 and beta_true_gen
  save(result_matrix_adjmlecoverage, file=paste0("result_matrix_adjmlecoverage_kappa1_",kappa1, "_delta_",delta,"_beta_true_gen_",beta_true_gen,".RData"))
  
}

for(kappa1 in c(0.5,1,1.5,2)){
  for(delta in c(4)){
    coverage_adjMLE(kappa1,delta,beta_true_gen="t3",N_rep=50,M_div_p=20,p_seq=c(100 , 400 ,1600))
  }
}




coverage_MLE<-function(kappa1,delta,beta_true_gen="t3",N_rep=100,M_div_p=20,p_seq=c(100 , 400 ,1600)){
  result_matrix_mlecoverage=matrix(0,nc=length(p_seq)-1,nr=N_rep)
  
  repeat_function_for_eta_square<-function(rep_i){
    cover_result_all_p=rep(0,length(p_seq)-1)
    set.seed(rep_i+100) # for observed data beta
    beta_true_large=generate_beta(p_seq[3],kappa1,beta_true_gen)
    set.seed(rep_i+200) # for observed data X
    X_large=matrix(rnorm(p_seq[3]^2*delta , mean = 0, sd = 1), nc = p_seq[3])
    
    for(pindex in 1:(length(p_seq)-1)){
      set.seed(rep_i+1000*pindex)
      p=p_seq[pindex]
      n=p*delta
      M=M_div_p*p
      beta_true=beta_true_large[1:p]
      
      X = X_large[1:n,1:p]/sqrt(p)
      Y = rbinom(n, 1, 1 / (1 + exp(-X %*% beta_true)))
      have_warning=FALSE
      tryCatch({
        fit_glm=glm(Y~X-1,family = binomial(link = "logit"))
      }, warning = function(w) {
        have_warning<<-TRUE
      })
      if(have_warning==FALSE){
        est=summary(fit_glm)$coefficients[,1]
        sdd=summary(fit_glm)$coefficients[,2]
        upper_CI=as.numeric(est+qnorm(0.975)*sdd)
        lower_CI=as.numeric(est-qnorm(0.975)*sdd)
        cover_result_all_p[pindex]=mean((beta_true>lower_CI)&(beta_true<upper_CI))
      }else{
        cover_result_all_p[pindex]=999
      }
    }
    return(cover_result_all_p)
  }
  
  parallel_eta_list=mclapply(1:N_rep,repeat_function_for_eta_square,mc.cores=numCores)
  for(rep_i in 1:N_rep){
    results_etasquare <- parallel_eta_list[[rep_i]]
    result_matrix_mlecoverage[rep_i,]=results_etasquare
  }
  
  
  # save file with kappa1 delta and tau_0 and beta_true_gen
  save(result_matrix_mlecoverage, file=paste0("result_matrix_mlecoverage_kappa1_",kappa1, "_delta_",delta,"_beta_true_gen_",beta_true_gen,".RData"))
  
}

for(kappa1 in c(0.5,1,1.5,2)){
  for(delta in c(4)){
    coverage_MLE(kappa1,delta,beta_true_gen="t3",N_rep=50,M_div_p=20,p_seq=c(100 , 400 ,1600))
  }
}