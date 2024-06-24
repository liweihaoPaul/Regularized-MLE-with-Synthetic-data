
#### generate plot
library(ggplot2)
library(sigmoid)
tau_0_seq=seq(0.1,2,0.1)
one_plot_pred_deviance<-function(kappa1,tau_0_seq=seq(0.1,2,0.1),design_type="Gaussian",beta_true_gen="t3"){
  delta=2
  load(paste0("TheorySolution_non_infor_syn_data_delta_kappa1_",kappa1,"_delta_",delta,".RData"))
  TheorySolution2=TheorySolution
  Z_matrix<-matrix(rnorm(100000),ncol=2)
  pred_deviance_seq_each_tau_0=rep(0,length(tau_0_seq))
  for (tau_0_index in 1:length(tau_0_seq)){
    tau_0=tau_0_seq[tau_0_index]
    pred_deviance_seq_each_tau_0[tau_0_index]=-mean(sigmoid(kappa1*Z_matrix[,2])*log(sigmoid(abs(TheorySolution2[tau_0_index,1])*Z_matrix[,1]+TheorySolution2[tau_0_index,3]*kappa1*Z_matrix[,2]))
                                                    +(1-sigmoid(kappa1*Z_matrix[,2]))*log(1-sigmoid(abs(TheorySolution2[tau_0_index,1])*Z_matrix[,1]+TheorySolution2[tau_0_index,3]*kappa1*Z_matrix[,2])))
  }
  plot_data_set=data.frame(x=tau_0_seq,y_mean=rep(0,length(tau_0_seq)),y_upper=rep(0,length(tau_0_seq)),
                           y_lower=rep(0,length(tau_0_seq)),theory_limit=pred_deviance_seq_each_tau_0)
  # add empirical point
  m=20/delta
  load(paste0("empirical_predresult_matrix_non_infor_syn_data_delta_",
              delta,"_m_",m,"_kappa1_",
              kappa1,"_design_type_",
              design_type,"_beta_true_gen_",beta_true_gen,".RData"))
   
  plot_data_set$y_mean=-empirical_result_matrix[,1]
  plot_data_set$delta=delta
  dataset1=plot_data_set
  
  
  
  
  
  delta=4
  load(paste0("TheorySolution_non_infor_syn_data_delta_kappa1_",kappa1,"_delta_",delta,".RData"))
  TheorySolution2=TheorySolution
  pred_deviance_seq_each_tau_0=rep(0,length(tau_0_seq))
  for (tau_0_index in 1:length(tau_0_seq)){
    tau_0=tau_0_seq[tau_0_index]
    pred_deviance_seq_each_tau_0[tau_0_index]=-mean(sigmoid(kappa1*Z_matrix[,2])*log(sigmoid(abs(TheorySolution2[tau_0_index,1])*Z_matrix[,1]+TheorySolution2[tau_0_index,3]*kappa1*Z_matrix[,2]))
                                                    +(1-sigmoid(kappa1*Z_matrix[,2]))*log(1-sigmoid(abs(TheorySolution2[tau_0_index,1])*Z_matrix[,1]+TheorySolution2[tau_0_index,3]*kappa1*Z_matrix[,2])))
  }
  plot_data_set=data.frame(x=tau_0_seq,y_mean=rep(0,length(tau_0_seq)),y_upper=rep(0,length(tau_0_seq)),
                           y_lower=rep(0,length(tau_0_seq)),theory_limit=pred_deviance_seq_each_tau_0)
  
  # add empirical point
  m=20/delta
  load(paste0("empirical_predresult_matrix_non_infor_syn_data_delta_",
              delta,"_m_",m,"_kappa1_",
              kappa1,"_design_type_",
              design_type,"_beta_true_gen_",beta_true_gen,".RData"))
  plot_data_set$y_mean=-empirical_result_matrix[,1]
  plot_data_set$delta=delta
  dataset2=plot_data_set
  merged_data=rbind(dataset1,dataset2)
  figg=ggplot(merged_data, aes(x = x, group = delta)) +
    geom_line(aes(y = theory_limit, color = as.factor(delta))) +
    #geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = as.factor(delta)), alpha = 0.2, color = NA) +
    geom_point(aes(y = y_mean, color = as.factor(delta),shape= as.factor(delta)), size = 2) +
    scale_color_manual(values = c("red", "blue"), 
                       labels = c(expression(delta==2), expression(delta==4))) +
    scale_shape_manual(values = c(1,2),  # Assign specific shapes for each delta
                       labels = c(expression(delta == 2), expression(delta == 4))) +
    labs(x = expression(tau[0]), y = "Predictive Deviance", color = "" ) +
    ylim(0.5, 1) +
    theme(panel.background = element_blank(),  # Blank background
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          panel.border = element_rect(colour = "black", fill = NA) # Add the frame
          ,legend.position="none",plot.title = element_text(hjust = 0.5))+ggtitle(bquote(kappa[1] == .(kappa1))) 
  
  only_for_legend_plot=ggplot(merged_data, aes(x = x, group = delta)) +
    geom_line(aes(y = theory_limit, color = as.factor(delta))) +
    #geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = as.factor(delta)), alpha = 0.2, color = NA) +
    geom_point(aes(y = y_mean, color = as.factor(delta),shape= as.factor(delta)), size = 2) +
    scale_color_manual(values = c("red", "blue"),
                       labels = c(expression(delta==2), expression(delta==4))) +
    scale_shape_manual(values = c(1,2),  # Assign specific shapes for each delta
                       labels = c(expression(delta == 2), expression(delta == 4))) +
    labs(x = expression(tau[0]), y = "Cosine Similarity", color =  " ",shape= " " ) +
    theme(legend.position = "right")
 
  g_legend <- ggplotGrob(only_for_legend_plot)
  legend <- g_legend$grobs[[which(g_legend$layout$name == "guide-box")]]
  return(list(fig=figg,legend=legend))
  
}

one_plot_classerror<-function(kappa1,tau_0_seq=seq(0.1,2,0.1),design_type="Gaussian",beta_true_gen="t3"){
  delta=2
  load(paste0("TheorySolution_non_infor_syn_data_delta_kappa1_",kappa1,"_delta_",delta,".RData"))
  TheorySolution2=TheorySolution
  Z_matrix<-matrix(rnorm(100000),ncol=2)
  pred_deviance_seq_each_tau_0=rep(0,length(tau_0_seq))
  for (tau_0_index in 1:length(tau_0_seq)){
    tau_0=tau_0_seq[tau_0_index]
    pred_deviance_seq_each_tau_0[tau_0_index]=mean(sigmoid(kappa1*Z_matrix[,2])*(abs(TheorySolution2[tau_0_index,1])*Z_matrix[,1]+TheorySolution2[tau_0_index,3]*kappa1*Z_matrix[,2]<0)
                                                   +(1-sigmoid(kappa1*Z_matrix[,2]))*(abs(TheorySolution2[tau_0_index,1])*Z_matrix[,1]+TheorySolution2[tau_0_index,3]*kappa1*Z_matrix[,2]>0))
  }
  plot_data_set=data.frame(x=tau_0_seq,y_mean=rep(0,length(tau_0_seq)),y_upper=rep(0,length(tau_0_seq)),
                           y_lower=rep(0,length(tau_0_seq)),theory_limit=pred_deviance_seq_each_tau_0)
  # add empirical point
  m=20/delta
  load(paste0("empirical_predresult_matrix_non_infor_syn_data_delta_",
              delta,"_m_",m,"_kappa1_",
              kappa1,"_design_type_",
              design_type,"_beta_true_gen_",beta_true_gen,".RData"))
  plot_data_set$y_mean=empirical_result_matrix[,4]
  plot_data_set$delta=delta
  dataset1=plot_data_set
  
  
  
  
  delta=4
  load(paste0("TheorySolution_non_infor_syn_data_delta_kappa1_",kappa1,"_delta_",delta,".RData"))
  TheorySolution2=TheorySolution
  pred_deviance_seq_each_tau_0=rep(0,length(tau_0_seq))
  for (tau_0_index in 1:length(tau_0_seq)){
    tau_0=tau_0_seq[tau_0_index]
    pred_deviance_seq_each_tau_0[tau_0_index]=mean(sigmoid(kappa1*Z_matrix[,2])*(abs(TheorySolution2[tau_0_index,1])*Z_matrix[,1]+TheorySolution2[tau_0_index,3]*kappa1*Z_matrix[,2]<0)
                                                   +(1-sigmoid(kappa1*Z_matrix[,2]))*(abs(TheorySolution2[tau_0_index,1])*Z_matrix[,1]+TheorySolution2[tau_0_index,3]*kappa1*Z_matrix[,2]>0))
  }
  plot_data_set=data.frame(x=tau_0_seq,y_mean=rep(0,length(tau_0_seq)),y_upper=rep(0,length(tau_0_seq)),
                           y_lower=rep(0,length(tau_0_seq)),theory_limit=pred_deviance_seq_each_tau_0) 
  # add empirical point
  m=20/delta
  load(paste0("empirical_predresult_matrix_non_infor_syn_data_delta_",
              delta,"_m_",m,"_kappa1_",
              kappa1,"_design_type_",
              design_type,"_beta_true_gen_",beta_true_gen,".RData"))
  plot_data_set$y_mean=empirical_result_matrix[,4]
  plot_data_set$delta=delta
  dataset2=plot_data_set
  merged_data=rbind(dataset1,dataset2)
  figg=ggplot(merged_data, aes(x = x, group = delta)) +
    geom_line(aes(y = theory_limit, color = as.factor(delta))) +
    #geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = as.factor(delta)), alpha = 0.2, color = NA) +
    geom_point(aes(y = y_mean, color = as.factor(delta),shape= as.factor(delta)), size = 2) +
    scale_color_manual(values = c("red", "blue"), 
                       labels = c(expression(delta==2), expression(delta==4))) +
    scale_shape_manual(values = c(1,2),  # Assign specific shapes for each delta
                       labels = c(expression(delta == 2), expression(delta == 4))) +
    labs(x = expression(tau[0]), y = "Generalization Error", color = "" ) +
    ylim(0.3, 0.5) +
    theme(panel.background = element_blank(),  # Blank background
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          panel.border = element_rect(colour = "black", fill = NA) # Add the frame
          ,legend.position="none",plot.title = element_text(hjust = 0.5))+ggtitle(bquote(kappa[1] == .(kappa1))) 
  
  only_for_legend_plot=ggplot(merged_data, aes(x = x, group = delta)) +
    geom_line(aes(y = theory_limit, color = as.factor(delta))) +
    #geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = as.factor(delta)), alpha = 0.2, color = NA) +
    geom_point(aes(y = y_mean, color = as.factor(delta),shape= as.factor(delta)), size = 2) +
    scale_color_manual(values = c("red", "blue"),
                       labels = c(expression(delta==2), expression(delta==4))) +
    scale_shape_manual(values = c(1,2),  # Assign specific shapes for each delta
                       labels = c(expression(delta == 2), expression(delta == 4))) +
    labs(x = expression(tau[0]), y = "Cosine Similarity", color =  " ",shape= " " ) +
    theme(legend.position = "right",legend.text = element_text(size = 18),
          legend.key.size = unit(2.5, "lines"),panel.background = element_blank())
  
  g_legend <- ggplotGrob(only_for_legend_plot)
  legend <- g_legend$grobs[[which(g_legend$layout$name == "guide-box")]]
  return(list(fig=figg,legend=legend))
}

one_plot_classerror(kappa1=1)$fig
legend=one_plot_classerror(kappa1=1)$legend



library(gridExtra)
grid.arrange(one_plot_pred_deviance( kappa1 = 0.5)$fig, one_plot_pred_deviance( kappa1 = 1.5)$fig,
             one_plot_classerror( kappa1 = 0.5)$fig, one_plot_classerror( kappa1 = 1.5)$fig, nrow = 2,ncol=2,
             right=legend)

 



