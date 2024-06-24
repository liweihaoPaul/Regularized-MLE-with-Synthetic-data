
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

summary_one_setting<-function(kappa1,delta,kesee=0.9,kappa2=1,delta_s=10,
                              beta_true_gen="t3",N_rep=50,p_seq=c(100 , 400 ,1600)){
  load(paste0("parallel_betahatloocv_list_non_informative_syn_p400_kappa1_",
              kappa1,"_delta_",delta,".RData"))
  load(paste0("parallel_betahatestkappa1_mse_list_non_informative_syn_p400_kappa1_",
              kappa1,"_delta_",delta,".RData"))
  load(paste0("parallel_betahateoracle_mse_list_non_informative_syn_p400_kappa1_",
              kappa1,"_delta_",delta,".RData"))
  
  load(paste0("betahat_LOOCV_informative_syn_p400_kappa1_",kappa1,"_delta_",
              delta,"_kesee_",kesee,"_kappa2_",kappa2,"_delta_s_",delta_s,".RData"))
  load(paste0("betahat_estkappa1kappa2kesee_mse_curve_informative_syn_p400_kappa1_",kappa1,"_delta_",
              delta,"_kesee_",kesee,"_kappa2_",kappa2,"_delta_s_",delta_s,".RData"))
  load(paste0("parallel_betahatridgecv_list_non_informative_syn_p400_kappa1_",
              kappa1,"_delta_",delta,".RData"))
  load(paste0("betahat_oraclekesee_mse_curve_informative_syn_p400_kappa1_",
              kappa1,"_delta_",delta,"_kesee_",kesee,"_kappa2_",kappa2,"_delta_s_",delta_s,".RData"))
  mse_different_method=matrix(0,nrow=N_rep,ncol=7)
  for(rep_i in 1:N_rep){
    set.seed(rep_i+100) # for observed data beta
    beta_true_large=generate_beta(p_seq[3],kappa1,beta_true_gen)
    beta_true=beta_true_large[1:p_seq[2]]/sqrt(p_seq[2])
    mse_different_method[rep_i,]=c(norm(beta_true-parallel_betahatloocv_list[[rep_i]][-1],"2" )^2,
      norm(beta_true-parallel_betahatestkappa1_mse_list[[rep_i]][-1],"2" )^2,
      norm(beta_true-parallel_betahateoracle_mse_list[[rep_i]][-1],"2" )^2,
      norm(beta_true-parallel_betahatsourceloocv_list[[rep_i]][-1],"2" )^2,
      norm(beta_true-parallel_betahatsourceestkesee_mse_list[[rep_i]][-1],"2" )^2,
      norm(beta_true-parallel_betahatsourceoraclekesee_mse_list[[rep_i]][-1],"2" )^2,
      norm(beta_true-parallel_betahatridgecv_list[[rep_i]][-1],"2" )^2)
    
  }
  print(apply(mse_different_method,2,mean))
  print(apply(mse_different_method,2,sd))/sqrt(N_rep)
  #print(mse_different_method[,7])
  #print(sapply(1:N_rep,function(index){parallel_betahateoracle_mse_list[[index]][1]}))
  #print(mse_different_method[,1])
  #print(sapply(1:N_rep,function(index){parallel_betahatsourceoraclekesee_mse_list[[index]][1]}))
  boxplot(mse_different_method[,1:6],main=paste0("n/p=",delta),xlab="Method",ylab="Square Error") #main=paste0("kappa1=",kappa1,",delta=",delta,",p400"))
  # add legends for box plot, 7 box plot total from left to right
  # legend("bottom",inset=c(-0.3,0),legend=c("1 NoInfor:LOOCV","2 NoInfor:k1hat","3 NoInfor:oracle","4 Source:LOOCV",
  #                                            "5 Source:estAll","6 Source:oracle","7 Ridge: CV"),
  #        cex=0.8,title="Method")
}
par(mfrow=c(2,2))

for(kappa1 in c(1,2)){
  for(delta in c(2,4)){
    summary_one_setting(kappa1,delta)
  }
}


############## use ggplot to remake
library(tidyverse)
library(ggplot2)
one_ggplot<-function(kappa1,delta,ylimlow,ylimupp,kesee=0.9,kappa2=1,delta_s=10,
                     beta_true_gen="t3",N_rep=50,p_seq=c(100 , 400 ,1600)){
  load(paste0("parallel_betahatloocv_list_non_informative_syn_p400_kappa1_",
              kappa1,"_delta_",delta,".RData"))
  load(paste0("parallel_betahatestkappa1_mse_list_non_informative_syn_p400_kappa1_",
              kappa1,"_delta_",delta,".RData"))
  load(paste0("parallel_betahateoracle_mse_list_non_informative_syn_p400_kappa1_",
              kappa1,"_delta_",delta,".RData"))
  
  load(paste0("betahat_LOOCV_informative_syn_p400_kappa1_",kappa1,"_delta_",
              delta,"_kesee_",kesee,"_kappa2_",kappa2,"_delta_s_",delta_s,".RData"))
  load(paste0("betahat_estkappa1kappa2kesee_mse_curve_informative_syn_p400_kappa1_",kappa1,"_delta_",
              delta,"_kesee_",kesee,"_kappa2_",kappa2,"_delta_s_",delta_s,".RData"))
  load(paste0("parallel_betahatridgecv_list_non_informative_syn_p400_kappa1_",
              kappa1,"_delta_",delta,".RData"))
  load(paste0("betahat_oraclekesee_mse_curve_informative_syn_p400_kappa1_",
              kappa1,"_delta_",delta,"_kesee_",kesee,"_kappa2_",kappa2,"_delta_s_",delta_s,".RData"))
  mse_different_method=matrix(0,nrow=N_rep,ncol=7)
  for(rep_i in 1:N_rep){
    set.seed(rep_i+100) # for observed data beta
    beta_true_large=generate_beta(p_seq[3],kappa1,beta_true_gen)
    beta_true=beta_true_large[1:p_seq[2]]/sqrt(p_seq[2])
    mse_different_method[rep_i,]=c(norm(beta_true-parallel_betahatloocv_list[[rep_i]][-1],"2" )^2,
                                   norm(beta_true-parallel_betahatestkappa1_mse_list[[rep_i]][-1],"2" )^2,
                                   norm(beta_true-parallel_betahateoracle_mse_list[[rep_i]][-1],"2" )^2,
                                   norm(beta_true-parallel_betahatsourceloocv_list[[rep_i]][-1],"2" )^2,
                                   norm(beta_true-parallel_betahatsourceestkesee_mse_list[[rep_i]][-1],"2" )^2,
                                   norm(beta_true-parallel_betahatsourceoraclekesee_mse_list[[rep_i]][-1],"2" )^2,
                                   norm(beta_true-parallel_betahatridgecv_list[[rep_i]][-1],"2" )^2)
    
  }
  data_for_plot=mse_different_method[,1:6]
  data_for_plot=data.frame(data_for_plot)
  data_for_plot=data_for_plot %>% gather(key="Method",value="Square Error")
  
  ggplot(data_for_plot, aes(x=Method, y=`Square Error`)) +
    geom_boxplot() +
    labs(
      title = bquote(n == .(delta) * "p, " ~ kappa[1] == .(kappa1)),
      x = "Method",
      y = "Square Error"
    ) +
    scale_x_discrete( labels = c("1", "2", "3", "4", "5", "6")) + 
    ylim(c(ylimlow, ylimupp)) +  # Setting Y limits after data layers but before themes
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
  
  
}
p1=one_ggplot(1,2,0,2)
p2=one_ggplot(1,4,0,1.3)
p3=one_ggplot(2,2,0,6)
p4=one_ggplot(2,4,0,3.5)

library(gridExtra)
grid.arrange(p1,p2,p3,p4,ncol=2)






