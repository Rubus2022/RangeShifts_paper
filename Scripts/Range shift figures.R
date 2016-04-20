library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(igraph)
library(dplyr)
library(betalink)
library(tidyr)
library(viridis)
library(vegan)

source("./Functions/plotting_functions.R")

options(scipen=999)

#Figure 1####
load("./Workspace/Range_shift_heatplots.RData")
nprey<-40
npred1<-24
npred2<-16

par(mfrow=c(2,3))
local_net_plot(Com_inits = Com_list$FW1[,1,],Com_final = Com_list$FW1[,1,55],Ints = Int_list$B3,trophic = T)
local_net_plot(Com_inits = X2[,2000,],Com_final = X2[,7000,105],Ints = B3,trophic = T)
local_net_plot(Com_inits = X3[,2000,],Com_final = X3[,7000,105],Ints = B3,trophic = T)




#Figure 2####
Net_turn_means<-Net_shift.df%>%
  group_by(Part,Scale,Dispersal,Community)%>%
  summarise_each(funs(mean(.,na.rm=T),sd(.,na.rm=T)))


ggplot(Net_turn_means,aes(x=Dispersal,y=WN_mean,color=Community, fill=Community))+
  facet_grid(Part~Scale,scales="free")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=WN_mean-WN_sd,ymax=WN_mean+WN_sd),alpha=0.2, color=NA)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_bw(base_size = 15)+
  removeGrid()+
  ylab("Network dissimilarity")+
  xlab("Dispersal")
ggsave("./Figures/Network turnover.pdf",width = 8,height = 8)


Net_ind_means<-Net_ind.df%>%
  group_by(Scale,Dispersal,Community)%>%
  summarise_each(funs(mean(.,na.rm=T),sd(.,na.rm=T)))


ggplot(Net_ind_means,aes(x=Dispersal,y=LD_mean,color=Community, fill=Community))+
  facet_grid(.~Scale,scales="free")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=LD_mean-LD_sd,ymax=LD_mean+LD_sd),alpha=0.2, color=NA)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_bw(base_size = 15)+
  removeGrid()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1))+
  ylab("Link density (similarity)")+
  xlab("Dispersal")
ggsave("./Figures/Link density compare.pdf",width = 12,height = 5)

ggplot(filter(Net_ind_means,Community=="Food web"),aes(x=Dispersal,y=Trophic_levels_mean))+
  facet_grid(.~Scale,scales="free")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=Trophic_levels_mean-Trophic_levels_sd,ymax=Trophic_levels_mean+Trophic_levels_sd),alpha=0.2, color=NA)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_bw(base_size = 15)+
  removeGrid()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1))+
  ylab("# of trophic levels (similarity)")+
  xlab("Dispersal")
ggsave("./Figures/Trophic levels compare.pdf",width = 12,height = 5)


Trophic_shift_means<-Net_shift_troph.df%>%
  group_by(Part,Scale,Dispersal,Trophic)%>%
  summarise_each(funs(mean(.,na.rm=T),sd(.,na.rm=T)))

ggplot(Trophic_shift_means,aes(x=Dispersal,y=WN_mean,color=Trophic,fill=Trophic))+
  facet_grid(Part~Scale, scales="free")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=WN_mean-WN_sd,ymax=WN_mean+WN_sd),alpha=0.2,color=NA)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_bw(base_size = 15)+
  removeGrid()+
  ylab("Network dissimilarity")+
  xlab("Dispersal")
ggsave("./Figures/Trophic turnover.pdf",width = 8,height = 8)


Shift_means<-Shift.df%>%
  group_by(Interactions,Dispersal)%>%
  summarise_each(funs(mean(.,na.rm=T),sd(.,na.rm=T)))%>%
  mutate(Trophic="Community")%>%
  ungroup()%>%
  mutate(Interactions=revalue(Interactions,c("Carnivores"="Predators")))%>%
  mutate(Trophic=replace(Trophic,Interactions=="Plants","Food web"))%>%
  mutate(Trophic=replace(Trophic,Interactions=="Herbivores","Food web"))%>%
  mutate(Trophic=replace(Trophic,Interactions=="Predators","Food web"))


ggplot(Shift_means,aes(x=Dispersal,y=Shift_sd_mean,color=Interactions,fill=Interactions))+
  #facet_grid(Part~Scale, scales="free")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=Shift_sd_mean-Shift_sd_sd,ymax=Shift_sd_mean+Shift_sd_sd),alpha=0.2,color=NA)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  facet_wrap(~Trophic)+
  theme_bw(base_size = 15)+
  removeGrid()+
  ylab("Range shift variation")+
  xlab("Dispersal")
ggsave("./Figures/Range shift variation.pdf",width = 8,height = 8)
