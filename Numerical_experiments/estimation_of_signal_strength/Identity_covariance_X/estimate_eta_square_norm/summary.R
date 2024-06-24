#summary of the estimation


return_kappa1_given_eta_square<-function(pred_eta_square,Dict){
  closed_index=which.min(abs(Dict[,2]-pred_eta_square))
  return(Dict[closed_index,1])
}


summary_one_setting<-function(delta,tau_0,beta_true_gen,kappa1){
  load(paste0("correspondence_eta_kappa1_non_infor_syn_data_tau_0_0.25_delta_",delta,".RData"))
  load(paste0("result_matrix_eta_square_kappa1_",kappa1, "_delta_",delta,"_tau_0_",tau_0,"_beta_true_gen_",beta_true_gen,".RData"))
  result_matrix_eta_square=result_matrix_eta_square     #[-c(11,49,63,84),]
  est_kappa1_matrix=result_matrix_eta_square
  for(i in 1:nrow(result_matrix_eta_square)){
    for(j in 1:ncol(result_matrix_eta_square)){
      est_kappa1_matrix[i,j]=return_kappa1_given_eta_square(result_matrix_eta_square[i,j],correspondence_eta_kappa1)
    }
  }
  error_relative_matrix=abs(est_kappa1_matrix-kappa1)#/kappa1
  # print(rbind(apply(error_relative_matrix,2,median),
  # apply(error_relative_matrix,2,IQR)))
  print(round(rbind(apply(error_relative_matrix,2,mean),
              apply(error_relative_matrix,2,sd)),3))
}

summary_one_setting(delta=2,tau_0=1/4,beta_true_gen="t3",kappa1=0.5)
summary_one_setting(delta=2,tau_0=1/4,beta_true_gen="t3",kappa1=1)
summary_one_setting(delta=2,tau_0=1/4,beta_true_gen="t3",kappa1=1.5)
summary_one_setting(delta=2,tau_0=1/4,beta_true_gen="t3",kappa1=2)

summary_one_setting(delta=4,tau_0=1/4,beta_true_gen="t3",kappa1=0.5)
summary_one_setting(delta=4,tau_0=1/4,beta_true_gen="t3",kappa1=1)
summary_one_setting(delta=4,tau_0=1/4,beta_true_gen="t3",kappa1=1.5)
summary_one_setting(delta=4,tau_0=1/4,beta_true_gen="t3",kappa1=2)


kappa1_one_setting<-function(delta,tau_0,beta_true_gen,kappa1){
  load(paste0("correspondence_eta_kappa1_non_infor_syn_data_tau_0_0.25_delta_",delta,".RData"))
  load(paste0("result_matrix_eta_square_kappa1_",kappa1, "_delta_",delta,"_tau_0_",tau_0,"_beta_true_gen_",beta_true_gen,".RData"))
  result_matrix_eta_square=result_matrix_eta_square     #[-c(11,49,63,84),]
  est_kappa1_matrix=result_matrix_eta_square
  for(i in 1:nrow(result_matrix_eta_square)){
    for(j in 1:ncol(result_matrix_eta_square)){
      est_kappa1_matrix[i,j]=return_kappa1_given_eta_square(result_matrix_eta_square[i,j],correspondence_eta_kappa1)
    }
  }
  save(est_kappa1_matrix,file=paste0("result_matrix_kappa1_",kappa1,
                                     "_delta_",delta,"_tau_0_",
                                     tau_0,"_beta_true_gen_",beta_true_gen,".RData"))
}

for(delta in c(2,4)){
  for(kappa1 in c(0.5,1,1.5,2)){
    kappa1_one_setting(delta=delta,tau_0=1/4,beta_true_gen="t3",kappa1=kappa1)
  }
}



