# est risk curve with kappa1hat,kappa2hat, keseehat

source("system_eq_solver_infor_source_data.R")
library(glmnet)
library(parallel)
numCores <- detectCores()
print(numCores)
tau_0_seq=seq(0.1,2,0.1)
numCores=min(50,detectCores())
estimated_solution_p400<-function(kappa1,delta,kappa2,delta_s,kesee,tau_0=1/4,beta_true_gen="t3",N_rep=50,M_div_p=20,p_seq=c(100 , 400 ,1600)){
  # load estimate for kappa1 and kappa2
  load(paste0("result_matrix_kappa2_kappa1_is_",kappa1,"_delta_",delta,"_tau_0_",
              tau_0,"_kesee_",kesee,"_kappa2_",
              kappa2,"_delta_s_",delta_s,"_beta_true_gen_",
              beta_true_gen,".RData")) # est_kappa2_matrix
  kappa2_hat_p400=est_kappa2_matrix[,2]
  load(paste0("result_matrix_kappa1_",kappa1,
              "_delta_",delta,"_tau_0_",
              tau_0,"_beta_true_gen_",beta_true_gen,".RData")) #est_kappa1_matrix
  kappa1_hat_p400=est_kappa1_matrix[,2]
  # load result for estimate keseehat
  load(paste0("est_matrix_kesee_hat_kappa1_",
              kappa1,"_delta_",delta,"_tau_0_",
              tau_0,"_kesee_",kesee,"_kappa2_",
              kappa2,"_delta_s_",delta_s,"_beta_true_gen_",
              beta_true_gen,".RData"))
  kesee_hat_p400=result_matrix_xi_hat[,2]
  
  
  repeat_function<-function(indexx){
    kappa1_hat=kappa1_hat_p400[indexx]
    kappa2_hat=kappa2_hat_p400[indexx]
    kesee_hat=kesee_hat_p400[indexx]
    # store in a matrix
    solution_system_matrix=matrix(0,nrow=length(tau_0_seq),ncol=4)
    for(i in 1:length(tau_0_seq)){
      tau_0=tau_0_seq[i]
      m=delta_s/delta
      solution_system_matrix[i,]=fixed_point_infor_source_data(delta,m,kappa1_hat,kappa2_hat,kesee_hat,tau_0)
    }
    return(solution_system_matrix)
  }
  result_list_solution_kappa12keseehat_p400=mclapply(1:length(kappa1_hat_p400),repeat_function,mc.cores=numCores)
  save(result_list_solution_kappa12keseehat_p400,file=paste0("result_list_solution_kappa12keseehat_p400_kappa1_",
                                                             kappa1,"_delta_",delta,"_tau_0_",
                                                             tau_0,"_kesee_",kesee,"_kappa2_",
                                                             kappa2,"_delta_s_",delta_s,"_beta_true_gen_",
                                                             beta_true_gen,".RData"))
  
}

estimated_solution_p400(kappa1=1,delta=4,kappa2=1,delta_s=10,kesee=0.9)
estimated_solution_p400(kappa1=1,delta=2,kappa2=1,delta_s=10,kesee=0.9)
estimated_solution_p400(kappa1=2,delta=4,kappa2=1,delta_s=10,kesee=0.9)
estimated_solution_p400(kappa1=2,delta=2,kappa2=1,delta_s=10,kesee=0.9)

