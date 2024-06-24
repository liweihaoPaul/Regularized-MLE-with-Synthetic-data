rm(list = ls())
library(MASS)
# library(mvnfast)
# library(knockoff)
library(glmhd)
library(parallel)
num_cores <- detectCores()
### Set working directory to 'figure2' folder, e.g., 
# setwd("~/code/simulation/figure2")

### source code
# source('utils.R')
# source('DS.R')
# source('MDS.R')
# source('BHq.R')
# source('GM.R')
# source('ABHq.R')
r_files <- list.files(pattern = "\\.R$", full.names = TRUE)
r_files<-r_files[-which(r_files %in% c("./vary_cor.R","qqqq"))]

for (file in r_files) {
  source(file)
}


### replicate index
#replicate <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#set.seed(replicate)

f <- function(x){
  exp(x)/(1 + exp(x))
}

p = 60
n = 500
s = 30
q = 0.1



for(rho in c(0,0.1,0.2,0.3,0.4)){
   repeat_function<-function(rep_i){
     set.seed(rep_i)
     signal_strength = 6.5
     Sigma = matrix(0, nrow = p, ncol = p)
     for(i in 1:p){
       for(j in 1:p){
         Sigma[i,j] = rho^abs(i-j)
       }
     }
     X = mvrnorm(n, mu = rep(0, p), Sigma = Sigma)*1/sqrt(n)
     signal_index = sample(1:p, size = s, replace = F)
     beta = numeric(p)
     beta[signal_index] = sample(c(-1, 1), s, replace = T)*signal_strength
     mu = X %*% beta
     y = rbinom(n, size = 1, prob = f(mu)) 
     
     # DS_time = system.time(DS_result <- DS(X, y))[3]
     # DS_result <- fdp_power(DS_result$select_index)
     # DS_fdp <- DS_result$fdp
     # DS_power <- DS_result$power
     
     MDS_time = system.time(MDS_result <- MDS(X, y, num_split = 50))[3]
     MDS_result <- fdp_power(MDS_result$select_index, signal_index)
     MDS_fdp <- MDS_result$fdp
     MDS_power <- MDS_result$power
     
     BHq_time = system.time(BHq_result <- BHq(X, y))[3]
     BHq_result <- fdp_power(BHq_result$select_index, signal_index)
     BHq_fdp <- BHq_result$fdp
     BHq_power <- BHq_result$power
     
     GM_time = system.time(GM_result <- GM(X, y))[3]
     GM_result <- fdp_power(GM_result$select_index, signal_index)
     GM_fdp <- GM_result$fdp
     GM_power <- GM_result$power
     print("hhhh")
     ABHq_time = system.time(ABHq_result <- ABHq(X, y))[3]
     ABHq_result <- fdp_power(ABHq_result$select_index, signal_index)
     ABHq_fdp <- ABHq_result$fdp
     ABHq_power <- ABHq_result$power
     
     ### save data
     data_save <- list(#DS_fdp   = DS_fdp,   DS_power = DS_power, DS_time = DS_time,
       repindex = rep_i, rho = rho, p = p, n = n, s = s, q = q,
       MDS_fdp  = MDS_fdp,  MDS_power = MDS_power, MDS_time = MDS_time,
       BHq_fdp  = BHq_fdp,  BHq_power = BHq_power, BHq_time = BHq_time,
       GM_fdp   = GM_fdp,   GM_power = GM_power, GM_time = GM_time,
       ABHq_fdp = ABHq_fdp, ABHq_power = ABHq_power, ABHq_time = ABHq_time)
     return(data_save)
   }
   data_save_list=mclapply(1:50, repeat_function, mc.cores = num_cores)
   save(data_save_list, file = paste0('vary_cor_',rho,'.RData'))
}





# summary
for(rho in c(0,0.1,0.2,0.3,0.4)){
  load(paste0('vary_cor_',rho,'.RData'))
  result_matrix=matrix(unlist(data_save_list),nc=18,byrow = T)
  fdr_col_index=c(7,10,13,16)
  power_col_index=c(8,11,14,17)
  print(rho)
  print( apply(result_matrix[,fdr_col_index],2,mean))
  print( apply(result_matrix[,power_col_index],2,mean))
}







