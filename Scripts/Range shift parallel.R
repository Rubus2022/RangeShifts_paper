library(vegan)
library(dplyr)
library(tidyr)
library(betalink)
library(igraph)
library(NetIndices)
library(RInSp)
#library(doParallel)
#library(foreach)
library(parallel)

#source("./Functions/process_nets.R")
source("./Functions/rShift_sim - 2.1.R")


#run simulation function in parallel
reps<-2
Sim_data_parallel<-mclapply(X = 1:reps, FUN= rShift_sim, mc.cores = 2)

# Sim_data_parallel<-foreach(r = 1:reps,.packages=c("igraph","dplyr","tidyr","RInSp","vegan","betalink","NetIndices")) %dopar% rShift_sim()
# stopCluster()
for(r in 1:reps){
  Sim_data<-Sim_data_parallel[[r]]
  Net_shift.df_temp<-Sim_data[[1]]
  Net_shift_troph.df_temp<-Sim_data[[2]]
  Net_ind.df_temp<-Sim_data[[3]]
  Net_ind_troph.df_temp<-Sim_data[[4]]
  Shift.df_temp<-Sim_data[[5]]
  if(r==1){
    Net_shift.df<-Net_shift.df_temp
    Net_shift_troph.df<-Net_shift_troph.df_temp
    Net_ind.df<-Net_ind.df_temp
    Net_ind_troph.df<-Net_ind_troph.df_temp
    Shift.df<-Shift.df_temp
  } else {
    Net_shift.df<-rbind(Net_shift.df,Net_shift.df_temp)
    Net_shift_troph.df<-rbind(Net_shift_troph.df,Net_shift_troph.df_temp)
    Net_ind.df<-rbind(Net_ind.df,Net_ind.df_temp)
    Net_ind_troph.df<-rbind(Net_ind_troph.df,Net_ind_troph.df_temp)
    Shift.df<-rbind(Shift.df,Shift.df_temp)
  }
}

save(Net_shift.df,Net_shift_troph.df,Net_ind.df,Net_ind_troph.df,Shift.df,file="Range_shift.RData")

