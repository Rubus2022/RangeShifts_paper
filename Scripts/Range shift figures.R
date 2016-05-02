library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(igraph)
library(betalink)
library(tidyr)
library(viridis)
library(vegan)
library(dplyr)
library(gridExtra)
library(cowplot)

source("./Functions/plotting_functions.R")

options(scipen=999)

load("./Workspace/Range_shift.RData")
Net_turn_means<-Net_shift.df%>%
  group_by(Part,Scale,Dispersal,Community)%>%
  summarise_each(funs(mean(.,na.rm=T),lower=quantile(.,probs=c(0.25),na.rm=T),upper=quantile(.,probs=c(0.75),na.rm=T)))

#Fig 1 bar charts####
pdf("./Figures/Fig. 1 barcharts.pdf", width = 10,height = 8)
par(mfrow=c(3,4),mar=c(1,1,1,1),oma=c(2,4,4,2))
hold<-filter(Net_turn_means,Scale=="Local network", Community=="Competitive",Dispersal==0.001)$WN_mean
hold[1]<-1-hold[1]
barplot(matrix(hold,ncol = 1),beside = F, col=c("dodgerblue3", "#ff7f00","grey"),horiz = T,add=F,width = 0.1, ylim=c(0,1),xaxt='n', ann=FALSE)

hold<-filter(Net_turn_means,Scale=="Local network", Community=="Competitive",Dispersal==0.01)$WN_mean
hold[1]<-1-hold[1]
barplot(matrix(hold,ncol = 1),beside = F, col=c("dodgerblue3", "#ff7f00","grey"),horiz = T,add=F,width = 0.1, ylim=c(0,1),xaxt='n', ann=FALSE)

hold<-filter(Net_turn_means,Scale=="Local network", Community=="Competitive",Dispersal==0.1)$WN_mean
hold[1]<-1-hold[1]
barplot(matrix(hold,ncol = 1),beside = F, col=c("dodgerblue3", "#ff7f00","grey"),horiz = T,add=F,width = 0.1, ylim=c(0,1),xaxt='n', ann=FALSE)

hold<-filter(Net_turn_means,Scale=="Local network", Community=="Competitive",Dispersal==0.5)$WN_mean
hold[1]<-1-hold[1]
barplot(matrix(hold,ncol = 1),beside = F, col=c("dodgerblue3", "#ff7f00","grey"),horiz = T,add=F,width = 0.1, ylim=c(0,1),xaxt='n', ann=FALSE)  

hold<-filter(Net_turn_means,Scale=="Local network", Community=="Mixed",Dispersal==0.001)$WN_mean
hold[1]<-1-hold[1]
barplot(matrix(hold,ncol = 1),beside = F, col=c("dodgerblue3", "#ff7f00","grey"),horiz = T,add=F,width = 0.1, ylim=c(0,1),xaxt='n', ann=FALSE)

hold<-filter(Net_turn_means,Scale=="Local network", Community=="Mixed",Dispersal==0.01)$WN_mean
hold[1]<-1-hold[1]
barplot(matrix(hold,ncol = 1),beside = F, col=c("dodgerblue3", "#ff7f00","grey"),horiz = T,add=F,width = 0.1, ylim=c(0,1),xaxt='n', ann=FALSE)

hold<-filter(Net_turn_means,Scale=="Local network", Community=="Mixed",Dispersal==0.1)$WN_mean
hold[1]<-1-hold[1]
barplot(matrix(hold,ncol = 1),beside = F, col=c("dodgerblue3", "#ff7f00","grey"),horiz = T,add=F,width = 0.1, ylim=c(0,1),xaxt='n', ann=FALSE)

hold<-filter(Net_turn_means,Scale=="Local network", Community=="Mixed",Dispersal==0.5)$WN_mean
hold[1]<-1-hold[1]
barplot(matrix(hold,ncol = 1),beside = F, col=c("dodgerblue3", "#ff7f00","grey"),horiz = T,add=F,width = 0.1, ylim=c(0,1),xaxt='n', ann=FALSE) 

hold<-filter(Net_turn_means,Scale=="Local network", Community=="Food web",Dispersal==0.001)$WN_mean
hold[1]<-1-hold[1]
barplot(matrix(hold,ncol = 1),beside = F, col=c("dodgerblue3", "#ff7f00","grey"),horiz = T,add=F,width = 0.1, ylim=c(0,1),xaxt='n', ann=FALSE)

hold<-filter(Net_turn_means,Scale=="Local network", Community=="Food web",Dispersal==0.01)$WN_mean
hold[1]<-1-hold[1]
barplot(matrix(hold,ncol = 1),beside = F, col=c("dodgerblue3", "#ff7f00","grey"),horiz = T,add=F,width = 0.1, ylim=c(0,1),xaxt='n', ann=FALSE)

hold<-filter(Net_turn_means,Scale=="Local network", Community=="Food web",Dispersal==0.1)$WN_mean
hold[1]<-1-hold[1]
barplot(matrix(hold,ncol = 1),beside = F, col=c("dodgerblue3", "#ff7f00","grey"),horiz = T,add=F,width = 0.1, ylim=c(0,1),xaxt='n', ann=FALSE)

hold<-filter(Net_turn_means,Scale=="Local network", Community=="Food web",Dispersal==1)$WN_mean
hold[1]<-1-hold[1]
barplot(matrix(hold,ncol = 1),beside = F, col=c("dodgerblue3", "#ff7f00","grey"),horiz = T,add=F,width = 0.1, ylim=c(0,1),xaxt='n', ann=FALSE)  
dev.off()

#Figure 1####
load("./Workspace/Range_shift_heatplots.RData")
nprey<-40
npred1<-24
npred2<-16
n<-80

patchI<-105
patchM<-105
patchT<-105
adj_value<-0

#Fig 1####
pdf("./Figures/Figure 1.pdf",height=8,width = 10)
par(mfrow=c(3,4),mar=c(1,1,1,1),oma=c(2,4,4,2))
local_net_plot(Com_inits = Com_list$`Comp 1`[,,1],Com_final = Com_list$`Comp 1`[,patchI,101],Ints = Int_list$BI,trophic = F)
mtext(expression(bold("a")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`Comp 2`[,,1],Com_final = Com_list$`Comp 2`[,patchI,101],Ints = Int_list$BI,trophic = F)
mtext(expression(bold("b")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`Comp 3`[,,1],Com_final = Com_list$`Comp 3`[,patchI,101],Ints = Int_list$BI,trophic = F)
mtext(expression(bold("c")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`Comp 4`[,,1],Com_final = Com_list$`Comp 4`[,patchI,101],Ints = Int_list$BI,trophic = F)
mtext(expression(bold("d")),side = 3,adj=adj_value)

local_net_plot(Com_inits = Com_list$`Mixed 1`[,,1],Com_final = Com_list$`Mixed 1`[,patchI,101],Ints = Int_list$BM,trophic = F)
mtext(expression(bold("e")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`Mixed 2`[,,1],Com_final = Com_list$`Mixed 2`[,patchI,101],Ints = Int_list$BM,trophic = F)
mtext(expression(bold("f")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`Mixed 3`[,,1],Com_final = Com_list$`Mixed 3`[,patchI,101],Ints = Int_list$BM,trophic = F)
mtext(expression(bold("g")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`Mixed 4`[,,1],Com_final = Com_list$`Mixed 4`[,patchI,101],Ints = Int_list$BM,trophic = F)
mtext(expression(bold("h")),side = 3,adj=adj_value)

local_net_plot(Com_inits = Com_list$`FW 1`[,,1],Com_final = Com_list$`FW 1`[,patchI,101],Ints = Int_list$B3,trophic = T)
mtext(expression(bold("i")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`FW 2`[,,1],Com_final = Com_list$`FW 2`[,patchI,101],Ints = Int_list$B3,trophic = T)
mtext(expression(bold("j")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`FW 3`[,,1],Com_final = Com_list$`FW 3`[,patchI,101],Ints = Int_list$B3,trophic = T)
mtext(expression(bold("k")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`FW 4`[,,1],Com_final = Com_list$`FW 4`[,patchI,101],Ints = Int_list$B3,trophic = T)
mtext(expression(bold("l")),side = 3,adj=adj_value)

mtext("Dispersal = 0.001",side=3,outer = T,adj = 0.04,padj = -1)
mtext("Dispersal = 0.01",side=3,outer = T,adj = 0.35,padj = -1)
mtext("Dispersal = 0.1",side=3,outer = T,adj = 0.65,padj = -1)
mtext("Dispersal = 0.5",side=3,outer = T,adj = 0.95,padj = -1)
mtext("Competition",side=2,outer = T,adj = 0.9,padj = -1)
mtext("Mixed",side=2,outer = T,adj = 0.5,padj = -1)
mtext("Food web",side=2,outer = T,adj = 0.15,padj = -1)
dev.off()

#Figure 2####
library(plyr)
Shift_means<-Shift.df%>%
  group_by(Interactions,Dispersal)%>%
  summarise_each(funs(mean(.,na.rm=T),lower=quantile(.,probs=c(0.25),na.rm=T),upper=quantile(.,probs=c(0.75),na.rm=T)))%>%
  mutate(Trophic="Community")%>%
  ungroup()%>%
  mutate(Interactions=revalue(Interactions,c("Carnivores"="Predators")))%>%
  mutate(Trophic=replace(Trophic,Interactions=="Plants","Food web"))%>%
  mutate(Trophic=replace(Trophic,Interactions=="Herbivores","Food web"))%>%
  mutate(Trophic=replace(Trophic,Interactions=="Predators","Food web"))

Shift_means$Interactions<-factor(Shift_means$Interactions,levels = c("No interactions","Competition","Mixed","Food web","Plants","Herbivores","Predators"),ordered = T)

Trophic_shift_means<-Net_shift_troph.df%>%
  group_by(Part,Scale,Dispersal,Trophic)%>%
  summarise_each(funs(mean(.,na.rm=T),lower=quantile(.,probs=c(0.25),na.rm=T),upper=quantile(.,probs=c(0.75),na.rm=T)))
  
Fig_2a<-ggplot(filter(Net_turn_means,Part=="All",Scale=="Local network"),aes(x=Dispersal,y=WN_mean,color=Community, fill=Community))+
  #facet_grid(Part~Scale,scales="free")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=WN_lower,ymax=WN_upper),alpha=0.2, color=NA)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set1",labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
  scale_fill_brewer(palette = "Set1",labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
  theme_bw(base_size = 12)+
  removeGrid()+
  ylab("Network dissimilarity")+
  xlab("Dispersal")+
  theme(legend.position=c(1,1),legend.justification=c(1,1))

Fig_2b<-ggplot(filter(Shift_means,Trophic=="Community"),aes(x=Dispersal,y=Shift_sd_mean,color=Interactions,fill=Interactions))+
  #facet_grid(Part~Scale, scales="free")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=Shift_sd_lower,ymax=Shift_sd_upper),alpha=0.2, color=NA)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set1",labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
  scale_fill_brewer(palette = "Set1",labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
  #facet_wrap(~Trophic)+
  theme_bw(base_size = 12)+
  removeGrid()+
  ylab("Range shift rate variation (sd)")+
  xlab("Dispersal")+
  theme(legend.position=c(1,1),legend.justification=c(1,1))

Fig_2c<-ggplot(filter(Trophic_shift_means,Part=="All",Scale=="Local network"),aes(x=Dispersal,y=WN_mean,color=Trophic,fill=Trophic))+
  #facet_grid(Part~Scale, scales="free")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=WN_lower,ymax=WN_upper),alpha=0.2, color=NA)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set2",name="")+
  scale_fill_brewer(palette = "Set2",name="")+
  theme_bw(base_size = 12)+
  removeGrid()+
  ylab("Trophic network dissimilarity")+
  xlab("Dispersal")+
  theme(legend.position=c(0,0),legend.justification=c(0,0))

Fig_2d<-ggplot(filter(Shift_means,Trophic=="Food web"),aes(x=Dispersal,y=Shift_sd_mean,color=Interactions,fill=Interactions))+
  #facet_grid(Part~Scale, scales="free")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=Shift_sd_lower,ymax=Shift_sd_upper),alpha=0.2, color=NA)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set2",name="")+
  scale_fill_brewer(palette = "Set2",name="")+
  #facet_wrap(~Trophic)+
  theme_bw(base_size = 12)+
  removeGrid()+
  ylab("Range shift rate variation (sd)")+
  xlab("Dispersal")+
  theme(legend.position=c(1,1),legend.justification=c(1,1))

plot_grid(Fig_2a,Fig_2b,Fig_2c,Fig_2d,labels=c("a","b","c","d"),ncol=2,nrow=2)
ggsave("./Figures/Figure 2.pdf",width = 10,height = 8,scale = 1)

#Figure 3####
Net_ind_long<-gather(Net_ind.df,key = Measure, value = Proportion,N:Trophic_levels)
Net_ind_long_means<-Net_ind_long%>%
  group_by(Dispersal,Community, Scale,Measure)%>%
  summarise_each(funs(mean(.,na.rm=T),lower=quantile(.,probs=c(0.25),na.rm=T),upper=quantile(.,probs=c(0.75),na.rm=T)))

Net_ind_long_means<-filter(Net_ind_long_means,Scale=="Local network")
Net_ind_long_means<-filter(Net_ind_long_means,Measure=="N" | Measure=="Ltot" | Measure=="LD"| Measure=="C"|Measure=="Cbar"| Measure=="Nestedness")

Net_ind_long_means<-Net_ind_long_means%>%
  mutate(Measure=replace(Measure,Measure=="N","# of species"),
         Measure=replace(Measure,Measure=="Ltot","# of links"),
         Measure=replace(Measure,Measure=="LD","Link density"),
         Measure=replace(Measure,Measure=="C","Connectance"),
         Measure=replace(Measure,Measure=="Cbar","Compartmentalization"))

Net_ind_long_means$Measure<-factor(Net_ind_long_means$Measure,levels=c("# of species","# of links","Link density","Connectance","Nestedness","Compartmentalization"),ordered = T)

ggplot(Net_ind_long_means,aes(x=Dispersal,y=Proportion_mean,color=Community, fill=Community))+
  facet_wrap(~Measure,scales="free_y")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=Proportion_lower,ymax=Proportion_upper),alpha=0.2, color=NA)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set1",name="")+
  scale_fill_brewer(palette = "Set1",name="")+
  theme_bw(base_size = 12)+
  removeGrid()+
  geom_hline(yintercept = 1, linetype=2)+
  scale_y_continuous(breaks = c(0,0.5,1,1.5,2,2.5))+
  ylab("Propotion of initial network")+
  xlab("Dispersal")+
  theme(legend.position=c(1,1),legend.justification=c(1,1.6))
ggsave("./Figures/Figure 3.pdf",width = 10,height = 6,scale = 1)

#Figure 4####
Net_inds$Disp_text<-paste("Dispersal = ", Net_inds$Dispersal, sep="")
Net_inds_stand<-Net_inds%>%
  ungroup()%>%
  group_by(Community)%>%
  mutate_each(funs(decostand(.,method="range",na.rm=T)),-Dispersal,-Disp_text,-patch,-time)

ggplot(Net_inds_stand,aes(y=patch,x=time,fill=LD))+
  geom_raster()+
  facet_grid(Community~Disp_text)+
  theme_bw(base_size = 12)+
  removeGrid()+
  scale_color_viridis(option = "D",name="")+
  scale_fill_viridis(option = "D",name="")+
  ylab("Patch")+
  xlab("Time")
ggsave(filename = "./Figures/Figure 4.png",width = 8,height = 8,dpi = 300)

#Figure S1####
Fig_S1a<-ggplot(filter(Net_turn_means,Part=="Gain",Scale=="Local network"),aes(x=Dispersal,y=WN_mean,color=Community, fill=Community))+
  #facet_grid(Part~Scale,scales="free")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=WN_mean-WN_sd,ymax=WN_mean+WN_sd),alpha=0.2, color=NA)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set1",labels=c("None", "Competition", "Mixed","Food web"),name="")+
  scale_fill_brewer(palette = "Set1",labels=c("None", "Competition", "Mixed","Food web"),name="")+
  theme_bw(base_size = 12)+
  removeGrid()+
  ylab("Network dissimilarity (gains)")+
  xlab("Dispersal")+
  theme(legend.position="none")

Fig_S1b<-ggplot(filter(Net_turn_means,Part=="Loss",Scale=="Local network"),aes(x=Dispersal,y=WN_mean,color=Community, fill=Community))+
  #facet_grid(Part~Scale,scales="free")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=WN_mean-WN_sd,ymax=WN_mean+WN_sd),alpha=0.2, color=NA)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set1",labels=c("None", "Competition", "Mixed","Food web"),name="")+
  scale_fill_brewer(palette = "Set1",labels=c("None", "Competition", "Mixed","Food web"),name="")+
  theme_bw(base_size = 12)+
  removeGrid()+
  ylab("Network dissimilarity (losses)")+
  xlab("Dispersal")+
  theme(legend.position=c(1,1),legend.justification=c(1,1))

Fig_S1c<-ggplot(filter(Trophic_shift_means,Part=="Gain",Scale=="Local network"),aes(x=Dispersal,y=WN_mean,color=Trophic,fill=Trophic))+
  #facet_grid(Part~Scale, scales="free")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=WN_mean-WN_sd,ymax=WN_mean+WN_sd),alpha=0.2,color=NA)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  theme_bw(base_size = 12)+
  removeGrid()+
  ylab("Network dissimilarity (gains)")+
  xlab("Dispersal")+
  theme(legend.position="none")

Fig_S1d<-ggplot(filter(Trophic_shift_means,Part=="Loss",Scale=="Local network"),aes(x=Dispersal,y=WN_mean,color=Trophic,fill=Trophic))+
  #facet_grid(Part~Scale, scales="free")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=WN_mean-WN_sd,ymax=WN_mean+WN_sd),alpha=0.2,color=NA)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set2",name="")+
  scale_fill_brewer(palette = "Set2",name="")+
  theme_bw(base_size = 12)+
  removeGrid()+
  ylab("Network dissimilarity (losses)")+
  xlab("Dispersal")+
  theme(legend.position=c(0,0),legend.justification=c(0,0))


plot_grid(Fig_S1a,Fig_S1b,Fig_S1c,Fig_S1d,labels=c("a)","b)","c)","d)"),ncol=2,nrow=2)
ggsave("./Figures/Figure S1.pdf",width = 10,height = 8,scale = 1)




