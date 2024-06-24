# native correlation between source and target data



library(glmnet)

return_alpha_given_kappa<-function(kappa1T,Theorysolution,kappa1_seq=seq(0.2,6,0.02)){
  closed_index=which.min(abs(kappa1_seq-kappa1T))
  return(Theorysolution[closed_index,3])
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

return_one_time_native_xi<-function(rep_i,kappa1,delta,tau_0,
                             kesee,kappa2,delta_s,beta_true_gen="t3",
                             N_rep=50,M_div_p=20,p_seq=c(100 , 400 ,1600)){
  three_nativexi_hat<-rep(0,3)
  betahat_obs_list<-list()
  set.seed(rep_i+100) # for observed data beta
  beta_true_large=generate_beta(p_seq[3],kappa1,beta_true_gen)
  set.seed(rep_i+200) # for observed data X
  X_large=matrix(rnorm(p_seq[3]^2*delta , mean = 0, sd = 1), nc = p_seq[3])
  # set.seed(rep_i+300) # for source data beta and X
  # betas_large=kesee*kappa2/kappa1*beta_true_large+lambda*generate_beta(p_seq[3],kappa1=1,beta_true_gen)
  # X_large=matrix(rnorm(p_seq[3]^2*delta_s , mean = 0, sd = 1), nc = p_seq[3]) # make sure generate different data
  set.seed(rep_i+400) # for synthetic data
  Xstar_large=matrix(rnorm(p_seq[3]^2*M_div_p , mean = 0, sd = 1), nc = p_seq[3])
  for(pindex in 1:length(p_seq)){
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
    #Ystar=Ystar_large[1:M]
    fit_cat_tau_0=glmnet(rbind(X,Xstar),c(Y,Ystar),weights =c(rep(1,n),rep(tau_0*n/M,M)) ,
                         family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
    betahat_tau0=as.numeric(fit_cat_tau_0$beta)
    betahat_obs_list[[pindex]]=betahat_tau0
  }
  betahat_source_list<-list()
  set.seed(rep_i+100) # for observed data beta
  beta_true_large=generate_beta(p_seq[3],kappa1,beta_true_gen)
  # set.seed(rep_i+200) # for observed data X
  # X_large=matrix(rnorm(p_seq[3]^2*delta , mean = 0, sd = 1), nc = p_seq[3])
  set.seed(rep_i+300) # for source data beta and X
  lambda=kappa2*(sqrt(1-kesee^2))
  betas_large=kesee*kappa2/kappa1*beta_true_large+lambda*generate_beta(p_seq[3],kappa1=1,beta_true_gen)
  X_large=matrix(rnorm(p_seq[3]^2*delta_s , mean = 0, sd = 1), nc = p_seq[3]) # make sure generate different data
  set.seed(rep_i+400) # for synthetic data
  Xstar_large=matrix(rnorm(p_seq[3]^2*M_div_p , mean = 0, sd = 1), nc = p_seq[3])
  for(pindex in 1:length(p_seq)){
    set.seed(rep_i+1000*pindex)
    p=p_seq[pindex]
    n=p*delta_s
    M=M_div_p*p
    betas=betas_large[1:p]/sqrt(p)
    
    X = X_large[1:n,1:p]
    Y = rbinom(n, 1, 1 / (1 + exp(-X %*% betas)))
    Xstar = Xstar_large[1:M,1:p]
    beta0 = rep(0, p)
    Ystar = rbinom(M, 1, 1 / (1 + exp(-Xstar %*% beta0)))
    fit_cat_tau_0=glmnet(rbind(X,Xstar),c(Y,Ystar),weights =c(rep(1,n),rep(tau_0*n/M,M)) ,
                         family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
    betahat_tau0=as.numeric(fit_cat_tau_0$beta) 
    betahat_source_list[[pindex]]=betahat_tau0
  }
  
  for(pindex in 1:length(p_seq)){
    
    three_nativexi_hat[pindex]=cor(betahat_obs_list[[pindex]],betahat_source_list[[pindex]])   #min(1,sum(betahat_obs_list[[pindex]]*betahat_source_list[[pindex]])/alpha1hat/alpha2hat/kappa1hat/kappa2hat)
  }
  
  return(three_nativexi_hat)
}

library(parallel)
numCores <- detectCores()

est_nativexiRelaerror_one_setting<-function(kappa1,delta,tau_0=1/4,kesee,kappa2,delta_s,beta_true_gen="t3",N_rep=50,M_div_p=20,p_seq=c(100 , 400 ,1600)){
  # use mclapply
  all_xihat_list=mclapply(1:N_rep, function(rep_i) {
    return_one_time_native_xi(rep_i,kappa1,delta,tau_0,kesee,kappa2,delta_s,beta_true_gen,N_rep,M_div_p,p_seq)
  }, mc.cores = numCores)
  result_matrix_xi_hat=do.call(rbind,all_xihat_list)
  error_relative_matrix=abs(result_matrix_xi_hat-kesee)/kesee
  save(error_relative_matrix,file=paste0("error_matrix_nativekesee_hat_kappa1_",
                                         kappa1,"_delta_",delta,"_tau_0_",
                                         tau_0,"_kesee_",kesee,"_kappa2_",
                                         kappa2,"_delta_s_",delta_s,"_beta_true_gen_",
                                         beta_true_gen,".RData"))
}


for(kappa1 in c(1,2)){
  for(delta in c(2,4)){
    for(kappa2 in c(1,2)){
      for(delta_s in c(4,10)){
        for(kesee in c(0.5,0.9)){
          est_nativexiRelaerror_one_setting(kappa1,delta,tau_0=1/4,kesee,kappa2,delta_s,beta_true_gen="t3",N_rep=50,M_div_p=20,p_seq=c(100 , 400 ,1600))
        }
      }
    }
  }
}


tau_0=1/4
beta_true_gen="t3"
for(kappa1 in c(1)){
  for(delta in c(2,4)){
    for(kappa2 in c(1)){
      for(delta_s in c(4,10)){
        for(kesee in c(0.9)){
          load(paste0("error_matrix_nativekesee_hat_kappa1_",
                 kappa1,"_delta_",delta,"_tau_0_",
                 tau_0,"_kesee_",kesee,"_kappa2_",
                 kappa2,"_delta_s_",delta_s,"_beta_true_gen_",
                 beta_true_gen,".RData"))
          print_matrix=matrix(0,nc=3,nr=2)
          error_relative_matrix=error_relative_matrix*kesee
          print_matrix[1,]=apply(error_relative_matrix,2,mean)
          print_matrix[2,]=apply(error_relative_matrix,2,sd)
          print_matrix=round(print_matrix,3)
          colnames(print_matrix)=c("p=100","p=400","p=1600")
          rownames(print_matrix)=c("mean","sd")
          # print current setting
          print(paste0("kappa1=",kappa1,"delta=",delta,"kappa2=",kappa2,"delta_s=",delta_s,"kesee=",kesee))
          print(print_matrix)
        }
      }
    }
  }
}



