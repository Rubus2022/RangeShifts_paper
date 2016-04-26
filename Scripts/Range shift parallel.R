library(vegan)
library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(betalink)
library(igraph)
library(viridis)
library(NetIndices)
library(RInSp)
library(doParallel)
library(foreach)

#source("./Functions/process_nets.R")
source("./Functions/rShift_sim.r")


#set up parallel####
cl<-makeCluster(detectCores())
registerDoParallel(cl)
getDoParWorkers()

#simulation code####
reps<-2

#run simulation function in parallel
Sim_data_parallel<-foreach(r = 1:reps,.packages=c("igraph","dplyr","tidyr","vegan","betalink","NetIndices")) %dopar% rShift_sim()
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

