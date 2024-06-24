
library(ggplot2)
library(reshape2)
library(cowplot)
library(gridExtra)
fdrcor_dataframe <- data.frame(cor=c(0,0.2,0.4,0.6,0.8), m1=c(0,0.2,0.4,0.6,0.8),
                            m2=c(0,0.2,0.4,0.6,0.8),m3=c(0,0.2,0.4,0.6,0.8),
                            m4=c(0,0.2,0.4,0.6,0.8),m5=c(0,0.2,0.4,0.6,0.8),m6=c(0,0.2,0.4,0.6,0.8),m7=c(0,0.2,0.4,0.6,0.8))
powercor_dataframe <- data.frame(cor=c(0,0.2,0.4,0.6,0.8), m1=c(0,0.2,0.4,0.6,0.8),
                               m2=c(0,0.2,0.4,0.6,0.8),m3=c(0,0.2,0.4,0.6,0.8),
                               m4=c(0,0.2,0.4,0.6,0.8),m5=c(0,0.2,0.4,0.6,0.8),m6=c(0,0.2,0.4,0.6,0.8),m7=c(0,0.2,0.4,0.6,0.8))
for(rho in c(0,0.2,0.4,0.6,0.8)){
  load(paste0('p500n3000MAPvary_cor_',rho,'.RData'))
  result_matrix=matrix(unlist(data_save_list),nc=15,byrow = T)
  fdr_col_index=c(7,10,13 )
  power_col_index=c(8,11,14 )
  fdrcor_dataframe[which(rho==c(0,0.2,0.4,0.6,0.8)),c(6,7,8)]=apply(result_matrix[,fdr_col_index],2,mean)
  powercor_dataframe[which(rho==c(0,0.2,0.4,0.6,0.8)),c(6,7,8)]=apply(result_matrix[,power_col_index],2,mean)
  load(paste0('p500n3000vary_cor_',rho,'.RData'))
  result_matrix=matrix(unlist(data_save_list),nc=18,byrow = T)
  fdr_col_index=c(7,10,13,16)
  power_col_index=c(8,11,14,17)
  fdrcor_dataframe[which(rho==c(0,0.2,0.4,0.6,0.8)),c(2,3,4,5)]=apply(result_matrix[,fdr_col_index],2,mean)
  powercor_dataframe[which(rho==c(0,0.2,0.4,0.6,0.8)),c(2,3,4,5)]=apply(result_matrix[,power_col_index],2,mean)
}
long_dfpower <- melt(powercor_dataframe, id.vars = 'cor', variable.name = 'Method', value.name = 'Error')
long_dfpower$Method <-rep(c('MDS(MLE)', 'BH(MLE)', 'GM(MLE)','ABH(MLE)', 'MDS(MAP)', 'ABH(MAP)', 'ABY(MAP)'),each=5)
long_dffdr <- melt(fdrcor_dataframe, id.vars = 'cor', variable.name = 'Method', value.name = 'Error')
long_dffdr$Method <- rep(c('MDS(MLE)', 'BH(MLE)', 'GM(MLE)','ABH(MLE)', 'MDS(MAP)', 'ABH(MAP)', 'ABY(MAP)'),each=5)

plotcorup=ggplot(long_dffdr[-c(11:15),], aes(x = cor, y = Error, color = Method,shape=Method)) +
  geom_point() +  # or geom_point() if you prefer points
  #scale_shape_manual(values=1:7) +
  labs(y = "FDR",x="") +
  theme_minimal()+ylim(0.,0.25)
plotcordown=ggplot(long_dfpower[-c(11:15),], aes(x = cor, y = Error, color = Method,shape=Method)) +
  geom_point() +  # or geom_point() if you prefer points
  #scale_shape_manual(values=1:7) +
  labs(x = "Toeplitz correlation",
       y = "Power") +
  theme_minimal()+ylim(0,1)


######### now for signal strength part
fdrsig_dataframe <- data.frame(Signal=c(8,9,10,11,12), m1=c(0,0.2,0.4,0.6,0.8),
                               m2=c(0,0.2,0.4,0.6,0.8),m3=c(0,0.2,0.4,0.6,0.8),
                               m4=c(0,0.2,0.4,0.6,0.8),m5=c(0,0.2,0.4,0.6,0.8),m6=c(0,0.2,0.4,0.6,0.8),m7=c(0,0.2,0.4,0.6,0.8))
powersig_dataframe <- data.frame(Signal=c(8,9,10,11,12), m1=c(0,0.2,0.4,0.6,0.8),
                                 m2=c(0,0.2,0.4,0.6,0.8),m3=c(0,0.2,0.4,0.6,0.8),
                                 m4=c(0,0.2,0.4,0.6,0.8),m5=c(0,0.2,0.4,0.6,0.8),m6=c(0,0.2,0.4,0.6,0.8),m7=c(0,0.2,0.4,0.6,0.8))

for(signal_strengthscalr in c(8,9,10,11,12)){
  load( paste0('p500n3000MAPvary_signal_',signal_strengthscalr,'.RData'))
  result_matrix=matrix(unlist(data_save_list),nc=15,byrow = T)
  fdr_col_index=c(7,10,13 )
  power_col_index=c(8,11,14 )
  fdrsig_dataframe[which(signal_strengthscalr== c(8,9,10,11,12)),c(6,7,8)]=apply(result_matrix[,fdr_col_index],2,mean)
  powersig_dataframe[which(signal_strengthscalr == c(8,9,10,11,12)),c(6,7,8)]=apply(result_matrix[,power_col_index],2,mean)
  load( paste0('p500n3000vary_signal_',signal_strengthscalr,'.RData'))
  result_matrix=matrix(unlist(data_save_list),nc=18,byrow = T)
  fdr_col_index=c(7,10,13,16)
  power_col_index=c(8,11,14,17)
  fdrsig_dataframe[which(signal_strengthscalr== c(8,9,10,11,12)),c(2,3,4,5)]=apply(result_matrix[,fdr_col_index],2,mean)
  powersig_dataframe[which(signal_strengthscalr== c(8,9,10,11,12)),c(2,3,4,5)]=apply(result_matrix[,power_col_index],2,mean)
}
long_dfpower <- melt(powersig_dataframe, id.vars = 'Signal', variable.name = 'Method', value.name = 'Error')
long_dfpower$Method <- rep(c('MDS(MLE)', 'BH(MLE)','GM(MLE)', 'ABH(MLE)','MDS(MAP)', 'ABH(MAP)', 'ABY(MAP)'),each=5)
long_dffdr <- melt(fdrsig_dataframe, id.vars = 'Signal', variable.name = 'Method', value.name = 'Error')
long_dffdr$Method <- rep(c('MDS(MLE)', 'BH(MLE)', 'GM(MLE)','ABH(MLE)', 'MDS(MAP)', 'ABH(MAP)', 'ABY(MAP)'),each=5)

plotsigup=ggplot(long_dffdr[-c(11:15),], aes(x = Signal, y = Error, color = Method,shape=Method)) +
  geom_point() +  # or geom_point() if you prefer points
  #scale_shape_manual(values=1:7) +
  labs(y = "FDR",x="") +
  theme_minimal()+ylim(0.0,0.25)
plotsigdown=ggplot(long_dfpower[-c(11:15),], aes(x = Signal, y = Error, color = Method,shape=Method)) +
  geom_point() +  # or geom_point() if you prefer points
  #scale_shape_manual(values=1:7) +
  labs(x = "Signal strength",
       y = "Power") +
  theme_minimal()+ylim(0.45,1)


grid.arrange(plotcorup, plotsigup, plotcordown,plotsigdown, nrow = 2,ncol=2)



