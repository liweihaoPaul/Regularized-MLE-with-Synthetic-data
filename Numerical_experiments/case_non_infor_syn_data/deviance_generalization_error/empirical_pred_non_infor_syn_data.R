

library(glmnet)
library(sigmoid)
library(parallel)
numCores <- detectCores()

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
generate_design_X<-function(n,p,design_type="Gaussian"){
  if(design_type=="Gaussian"){
    return(matrix(rnorm(p*n,mean=0,sd=sqrt(1/p)),nc=p))
  }
  if(design_type=="t4"){
    return(matrix(rt(p*n,df=4)/sqrt(p*2),nc=p))
  }
  if(design_type=="t3"){
    return(matrix(rt(p*n,df=3)/sqrt(p*3),nc=p))
  }
}

empirical_deviance_predicterror_non_infor_syn_data<-function(delta,m,kappa1,tau_0,
                                                  design_type="Gaussian",beta_true_gen="normal",p=250,N_test=500,NUM_REPEAT = 100){
  n=p*delta
  M=m*n
  tau=tau_0*n
  
  deviance_vector=rep(0, NUM_REPEAT)
  predicterror_vector=rep(0, NUM_REPEAT)
  for(rep_i in 1:NUM_REPEAT){
    set.seed(rep_i)
    beta_true=generate_beta(p,kappa1,beta_true_gen)
    X=generate_design_X(n,p,design_type)
    Xstar=matrix(rnorm(p*M,mean=0,sd=sqrt(1/p)),nc=p)
    prop1= 1/(1+exp(-X%*%beta_true))
    y1=rbinom(n,1,prop1)
    prop2= 1/(1+exp(rep(0,M)))
    ystar=rbinom(M,1,prop2)
    wt=c(rep(1,n),rep(tau/M,M))
    fit_lambda_seq<-glmnet(rbind(X,Xstar),c(y1,ystar),weights = wt,
                           family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
    beta_hat=as.numeric(fit_lambda_seq$beta )
    # reset seed to generate new design and response
    set.seed(rep_i+N_test*5)
    X_test=generate_design_X(N_test,p,design_type)
    prop1_test= 1/(1+exp(-X_test%*%beta_true))
    y1_test=rbinom(N_test,1,prop1_test)
    pred_deviance=mean(y1_test*log(sigmoid( X_test%*%beta_hat))+(1-y1_test)*log(1-sigmoid( X_test%*%beta_hat)))
    pred_y=as.numeric(sigmoid( X_test%*%beta_hat)>0.5)
    pred_error=mean(y1_test!=pred_y)
    deviance_vector[rep_i]=pred_deviance
    predicterror_vector[rep_i]=pred_error
  }
  # return average result
  return(c(mean(deviance_vector),mean(predicterror_vector))) 
}

tau_0_seq=seq(0.1,2,0.1)
for (kappa1 in c(0.5,1,1.5)){
  for (delta in c(2,4)){
    m=20/delta
    design_type="Gaussian"
    beta_true_gen="t3"
    empirical_result_matrix=matrix(0,nrow=length(tau_0_seq),ncol= 2)
    rep_function<-function(tau_0_index){
      
      tau_0=tau_0_seq[tau_0_index]
      empirical_result=empirical_deviance_predicterror_non_infor_syn_data(delta,m,kappa1,tau_0,
                                                               design_type,beta_true_gen)
      return(empirical_result)
    }
    
    empirical_pred_list=mclapply(1:length(tau_0_seq),rep_function,mc.cores=numCores)
    for(tau_0_index in 1:length(tau_0_seq)){
 
      empirical_result=empirical_pred_list[[tau_0_index]]
      empirical_result_matrix[tau_0_index,]=empirical_result
    }
    save(empirical_result_matrix,file=paste0("empirical_predresult_matrix_non_infor_syn_data_delta_",
                                             delta,"_m_",m,"_kappa1_",
                                             kappa1,"_design_type_",
                                             design_type,"_beta_true_gen_",beta_true_gen,".RData"))
  }
}





















