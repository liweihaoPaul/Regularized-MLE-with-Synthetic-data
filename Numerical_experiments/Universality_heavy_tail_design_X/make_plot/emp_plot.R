


generate_one_plot<-function(design_type="t3"){
  degree_f=c(2.5,3,4,5)[which(c("t25","t3","t4","t5")==design_type)]
  delta=4
  m=20/delta
  kappa1=1
  tau_0_seq=seq(0.1,2,0.1)
  load(paste0("system_root_non_infor_syn_data_delta_",delta,"_m_",m,"_kappa1_",kappa1,".RData"))
  beta_true_gen="t3"
  load(paste0("empirical_result_matrix_non_infor_syn_data_delta_",
              delta,"_m_",m,"_kappa1_",
              kappa1,"_design_type_",
              design_type,"_beta_true_gen_",beta_true_gen,".RData"))
  plot_data_set=data.frame(x=tau_0_seq,y_mean=rep(0,length(tau_0_seq)),y_upper=rep(0,length(tau_0_seq)),
                           y_lower=rep(0,length(tau_0_seq)),theory_limit=rep(0,length(tau_0_seq)))
  plot_data_set$theory_limit=system_root[,1]^2+(system_root[,3]-1)^2*kappa1^2
  plot_data_set$y_mean=empirical_result_matrix[,1]
   
  
  figg=ggplot(plot_data_set, aes(x = x, y = theory_limit)) +
    geom_line(aes(y = theory_limit)) +
    geom_point(aes(y = y_mean), size = 2) +
    labs(x = expression(tau[0]), y = "Square Error"  ) +
    ylim(0.6, 2.0) +
    theme(panel.background = element_blank(),  # Blank background
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          panel.border = element_rect(colour = "black", fill = NA) # Add the frame
          ,legend.position="none",plot.title = element_text(hjust = 0.5))+ggtitle(bquote(df == .(degree_f))) 
  return(figg)
 
}
generate_one_plot("t3")



library(gridExtra)
grid.arrange(generate_one_plot("t25"),generate_one_plot("t3"),generate_one_plot("t4"),generate_one_plot("t5"), nrow = 2,ncol=2)















