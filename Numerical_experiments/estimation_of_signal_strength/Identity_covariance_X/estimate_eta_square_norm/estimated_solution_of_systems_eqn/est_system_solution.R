# in this file we will look at the accuracy of estimated solution of system based on estimated kappa1

accurary_alpha_sigma_one_settings<-function(kappa1,delta,tau_0=1/4,beta_true_gen="t3"){
  kappa1_seq=seq(0.2,6,0.02)
  load(paste0("result_matrix_kappa1_",kappa1,
              "_delta_",delta,"_tau_0_",
              tau_0,"_beta_true_gen_",beta_true_gen,".RData"))
  load(paste0("TheorySolution_non_infor_syn_data_tau_0_0.25_delta_",delta,".RData"))
  true_solution=TheorySolution[which(kappa1_seq==kappa1),]
  print(round(true_solution[c(3,1)],3))
  est_solution_list_p_seq=list(p100=est_kappa1_matrix,p400=est_kappa1_matrix,p1600=est_kappa1_matrix)
  for(pindex in 1:3){
    for(i in 1:nrow(est_kappa1_matrix)){
      est_kappa1=est_kappa1_matrix[i,pindex]
      est_solution=TheorySolution[which(kappa1_seq==est_kappa1),]
      est_solution_list_p_seq[[pindex]][i,]=est_solution
    }
  }
  est_solution_error_list_p_seq=est_solution_list_p_seq
  for(pindex in 1:3){
    for(i in 1:nrow(est_kappa1_matrix)){
      est_solution_error_list_p_seq[[pindex]][i,]= abs(true_solution-est_solution_list_p_seq[[pindex]][i,])
    }
  }
  error_mean=c(apply(est_solution_error_list_p_seq[[1]],2,mean),
               apply(est_solution_error_list_p_seq[[2]],2,mean),
               apply(est_solution_error_list_p_seq[[3]],2,mean))
  error_sd=c(apply(est_solution_error_list_p_seq[[1]],2,sd),
               apply(est_solution_error_list_p_seq[[2]],2,sd),
               apply(est_solution_error_list_p_seq[[3]],2,sd))
  #return(round(rbind(error_mean,error_sd)[,-c(2,5,8)],3))
}

accurary_alpha_sigma_one_settings(kappa1 = 0.5,delta=2)
accurary_alpha_sigma_one_settings(kappa1 = 1,delta=2)
accurary_alpha_sigma_one_settings(kappa1 = 1.5,delta=2)
accurary_alpha_sigma_one_settings(kappa1 = 2,delta=2)

accurary_alpha_sigma_one_settings(kappa1 = 0.5,delta=4)
accurary_alpha_sigma_one_settings(kappa1 = 1,delta=4)
accurary_alpha_sigma_one_settings(kappa1 = 1.5,delta=4)
accurary_alpha_sigma_one_settings(kappa1 = 2,delta=4)
