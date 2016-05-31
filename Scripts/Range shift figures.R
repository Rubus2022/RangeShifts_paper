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

ColV<-c("grey20",brewer.pal(4,"Set1")[c(2,3,4)])

load("./Workspace/Range_shift.RData")
Net_turn_means<-Net_shift.df%>%
  group_by(Part,Scale,Dispersal,Community)%>%
  summarise_each(funs(mean(.,na.rm=T),lower=quantile(.,probs=c(0.25),na.rm=T),upper=quantile(.,probs=c(0.75),na.rm=T)))

Net_turn_means3<-Net_turn_means%>%
  filter(Community=="Food web",Dispersal>=0.05)

Net_turn_means_t<-Net_turn_means%>%
  mutate(WN_mean=replace(WN_mean,Community=="Food web" & Dispersal>0.05,NA))


#Fig 1 bar charts####
pdf("./Figures/Fig. 1 barcharts 2.pdf", width = 10,height = 8)
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

Fpatch<-120
adj_value<-0

#Fig 1####
pdf("./Figures/Figure 1.pdf",height=8,width = 10)
par(mfrow=c(3,4),mar=c(1,1,1,1),oma=c(2,4,4,2))
local_net_plot(Com_inits = Com_list$`Comp 1`[,,1],Com_final = Com_list$`Comp 1`[,Fpatch,101],Ints = Int_list$BI,trophic = F)
mtext(expression(bold("a")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`Comp 2`[,,1],Com_final = Com_list$`Comp 2`[,Fpatch,101],Ints = Int_list$BI,trophic = F)
mtext(expression(bold("b")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`Comp 3`[,,1],Com_final = Com_list$`Comp 3`[,Fpatch,101],Ints = Int_list$BI,trophic = F)
mtext(expression(bold("c")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`Comp 4`[,,1],Com_final = Com_list$`Comp 4`[,Fpatch,101],Ints = Int_list$BI,trophic = F)
mtext(expression(bold("d")),side = 3,adj=adj_value)

local_net_plot(Com_inits = Com_list$`Mixed 1`[,,1],Com_final = Com_list$`Mixed 1`[,Fpatch,101],Ints = Int_list$BM,trophic = F)
mtext(expression(bold("e")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`Mixed 2`[,,1],Com_final = Com_list$`Mixed 2`[,Fpatch,101],Ints = Int_list$BM,trophic = F)
mtext(expression(bold("f")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`Mixed 3`[,,1],Com_final = Com_list$`Mixed 3`[,Fpatch,101],Ints = Int_list$BM,trophic = F)
mtext(expression(bold("g")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`Mixed 4`[,,1],Com_final = Com_list$`Mixed 4`[,Fpatch,101],Ints = Int_list$BM,trophic = F)
mtext(expression(bold("h")),side = 3,adj=adj_value)

local_net_plot(Com_inits = Com_list$`FW 1`[,,1],Com_final = Com_list$`FW 1`[,Fpatch,101],Ints = Int_list$B3,trophic = T)
mtext(expression(bold("i")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`FW 2`[,,1],Com_final = Com_list$`FW 2`[,Fpatch,101],Ints = Int_list$B3,trophic = T)
mtext(expression(bold("j")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`FW 3`[,,1],Com_final = Com_list$`FW 3`[,Fpatch,101],Ints = Int_list$B3,trophic = T)
mtext(expression(bold("k")),side = 3,adj=adj_value)
local_net_plot(Com_inits = Com_list$`FW 4`[,,1],Com_final = Com_list$`FW 4`[,Fpatch,101],Ints = Int_list$B3,trophic = T)
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
Shift.df$Interactions<-factor(Shift.df$Interactions,levels=c("Predators",levels(Shift.df$Interactions)))

Shift_means<-Shift.df%>%
  group_by(Interactions,Dispersal)%>%
  summarise_each(funs(mean(.,na.rm=T),lower=quantile(.,probs=c(0.25),na.rm=T),upper=quantile(.,probs=c(0.75),na.rm=T)))%>%
  mutate(Trophic="Community")%>%
  ungroup()%>%
  mutate(Interactions=replace(Interactions,Interactions=="Carnivores","Predators"))%>%
  mutate(Trophic=replace(Trophic,Interactions=="Plants","Food web"))%>%
  mutate(Trophic=replace(Trophic,Interactions=="Herbivores","Food web"))%>%
  mutate(Trophic=replace(Trophic,Interactions=="Predators","Food web"))

Shift_means$Interactions<-factor(Shift_means$Interactions,levels = c("No interactions","Competition","Mixed","Food web","Plants","Herbivores","Predators"),ordered = T)

Shift_means3<-Shift_means%>%
  filter(Interactions=="Food web",Dispersal>=0.05)

Shift_means1<-Shift_means%>%
  mutate(Shift_sd_mean=replace(Shift_sd_mean,Interactions=="Food web" & Dispersal>0.05,NA))

Trophic_shift_means<-Net_shift_troph.df%>%
  group_by(Part,Scale,Dispersal,Trophic)%>%
  summarise_each(funs(mean(.,na.rm=T),lower=quantile(.,probs=c(0.25),na.rm=T),upper=quantile(.,probs=c(0.75),na.rm=T)))

Fig_2a<-ggplot(filter(Net_turn_means_t,Part=="All",Scale=="Local network"),aes(x=Dispersal,y=WN_mean,color=Community, fill=Community))+
  #facet_grid(Part~Scale,scales="free")+
  geom_ribbon(aes(ymin=WN_lower,ymax=WN_upper),alpha=0.2, color=NA)+
  geom_line(size=1.2)+
  geom_line(data=filter(Net_turn_means3,Part=="All",Scale=="Local network"),aes(x=Dispersal,y=WN_mean),size=1,linetype=3)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_manual(values = ColV,labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
  scale_fill_manual(values = ColV,labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
  theme_bw(base_size = 12)+
  removeGrid()+
  ylab("Network dissimilarity")+
  xlab("Dispersal")+
  theme(legend.position=c(1,1),legend.justification=c(1,1))

Fig_2b<-ggplot(filter(Shift_means1,Trophic=="Community"),aes(x=Dispersal,y=Shift_sd_mean,color=Interactions,fill=Interactions))+
  #facet_grid(Part~Scale, scales="free")+
  geom_ribbon(aes(ymin=Shift_sd_lower,ymax=Shift_sd_upper),alpha=0.2, color=NA)+
  geom_line(size=1.2)+
  geom_line(data=filter(Shift_means3,Trophic=="Community"),aes(x=Dispersal,y=Shift_sd_mean),size=1,linetype=3)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_manual(values = ColV,labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
  scale_fill_manual(values = ColV,labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
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
Net_ind.df<-Net_ind.df%>%
  mutate(Trophic_levels=ifelse(Community!="Food web",NA,Trophic_levels))

Net_ind_long<-gather(Net_ind.df,key = Measure, value = Proportion,N:Trophic_levels)

Net_ind_long_means<-Net_ind_long%>%
  group_by(Dispersal,Community, Scale,Measure)%>%
  summarise_each(funs(mean(.,na.rm=T),lower=quantile(.,probs=c(0.25),na.rm=T),upper=quantile(.,probs=c(0.75),na.rm=T)))

Net_ind_long_means<-filter(Net_ind_long_means,Scale=="Local network")
Net_ind_long_means<-filter(Net_ind_long_means,Measure=="N" | Measure=="Ltot" | Measure=="LD"| Measure=="Trophic_levels"|Measure=="Cbar"| Measure=="Nestedness")

Net_ind_long_means<-Net_ind_long_means%>%
  mutate(Measure=replace(Measure,Measure=="N","# of species"),
         Measure=replace(Measure,Measure=="Ltot","# of links"),
         Measure=replace(Measure,Measure=="LD","Link density"),
         Measure=replace(Measure,Measure=="Trophic_levels","Trophic levels"),
         Measure=replace(Measure,Measure=="Cbar","Compartmentalization"))


Net_ind_long_means$Measure<-factor(Net_ind_long_means$Measure,levels=c("# of species","# of links","Link density","Nestedness","Compartmentalization","Trophic levels"),ordered = T)


Net_ind_long_means_3<-Net_ind_long_means%>%
  filter(Community=="Food web",Dispersal>=0.05)

Net_ind_long_means<-Net_ind_long_means%>%
  mutate(Proportion_mean=replace(Proportion_mean,Community=="Food web" & Dispersal>0.05,NA))

ggplot(Net_ind_long_means,aes(x=Dispersal,y=Proportion_mean,color=Community, fill=Community))+
  facet_wrap(~Measure,scales="free_y")+
  geom_ribbon(aes(ymin=Proportion_lower,ymax=Proportion_upper),alpha=0.2, color=NA)+
  geom_line(size=1.2)+
  geom_line(data=Net_ind_long_means_3,aes(x=Dispersal,y=Proportion_mean),size=1,linetype=3)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_manual(values = ColV,labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
  scale_fill_manual(values = ColV,labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
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


init_patch_fun<-function(Com_final,Com_inits,Ints){
  nets_bin<-apply(cbind(Com_final,Com_inits),2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    Int_strength_cut<-quantile(Int_strength[Int_strength>0],0.5)#mean(Int_strength[Int_strength>0])
    Int_strength[Int_strength<Int_strength_cut]<-0
    return(c(Int_strength))
  })
  link_dis<-(as.matrix(vegdist(t(nets_bin),method = "jaccard",binary = T))[,1])[-1]
  min_dist<-which(link_dis==min(link_dis))
  return(min_dist)
}

contrast_data<-rbind(filter(Net_inds,
                            patch==max(init_patch_fun(Com_final = Com_list$`Comp 1`[,Fpatch,101],
                                                      Com_inits = Com_list$`Comp 1`[,,1],Ints = Int_list$BI)),
                            time==2000,Community=="Competition",Dispersal==0.001),
                     filter(Net_inds,
                            patch==max(init_patch_fun(Com_final = Com_list$`Comp 2`[,Fpatch,101],
                                                      Com_inits = Com_list$`Comp 2`[,,1],Ints = Int_list$BI)),
                            time==2000,Community=="Competition",Dispersal==0.01),
                     filter(Net_inds,
                            patch==max(init_patch_fun(Com_final = Com_list$`Comp 3`[,Fpatch,101],
                                                      Com_inits = Com_list$`Comp 3`[,,1],Ints = Int_list$BI)),
                            time==2000,Community=="Competition",Dispersal==0.1),
                     filter(Net_inds,
                            patch==max(init_patch_fun(Com_final = Com_list$`Comp 4`[,Fpatch,101],
                                                      Com_inits = Com_list$`Comp 4`[,,1],Ints = Int_list$BI)),
                            time==2000,Community=="Competition",Dispersal==0.5),
                     filter(Net_inds,
                            patch==max(init_patch_fun(Com_final = Com_list$`Mixed 1`[,Fpatch,101],
                                                      Com_inits = Com_list$`Mixed 1`[,,1],Ints = Int_list$BM)),
                            time==2000,Community=="Mixed interactions",Dispersal==0.001),
                     filter(Net_inds,
                            patch==max(init_patch_fun(Com_final = Com_list$`Mixed 2`[,Fpatch,101],
                                                      Com_inits = Com_list$`Mixed 2`[,,1],Ints = Int_list$BM)),
                            time==2000,Community=="Mixed interactions",Dispersal==0.01),
                     filter(Net_inds,
                            patch==max(init_patch_fun(Com_final = Com_list$`Mixed 3`[,Fpatch,101],
                                                      Com_inits = Com_list$`Mixed 3`[,,1],Ints = Int_list$BM)),
                            time==2000,Community=="Mixed interactions",Dispersal==0.1),
                     filter(Net_inds,
                            patch==max(init_patch_fun(Com_final = Com_list$`Mixed 4`[,Fpatch,101],
                                                      Com_inits = Com_list$`Mixed 4`[,,1],Ints = Int_list$BM)),
                            time==2000,Community=="Mixed interactions",Dispersal==0.5),
                     filter(Net_inds,
                            patch==max(init_patch_fun(Com_final = Com_list$`FW 1`[,Fpatch,101],
                                                      Com_inits = Com_list$`FW 1`[,,1],Ints = Int_list$B3)),
                            time==2000,Community=="Food web",Dispersal==0.001),
                     filter(Net_inds,
                            patch==max(init_patch_fun(Com_final = Com_list$`FW 2`[,Fpatch,101],
                                                      Com_inits = Com_list$`FW 2`[,,1],Ints = Int_list$B3)),
                            time==2000,Community=="Food web",Dispersal==0.01),
                     filter(Net_inds,
                            patch==max(init_patch_fun(Com_final = Com_list$`FW 3`[,Fpatch,101],
                                                      Com_inits = Com_list$`FW 3`[,,1],Ints = Int_list$B3)),
                            time==2000,Community=="Food web",Dispersal==0.1),
                     filter(Net_inds,
                            patch==max(init_patch_fun(Com_final = Com_list$`FW 4`[,Fpatch,101],
                                                      Com_inits = Com_list$`FW 4`[,,1],Ints = Int_list$B3)),
                            time==2000,Community=="Food web",Dispersal==0.5))


ggplot(Net_inds,aes(y=patch,x=time,fill=LD))+
  geom_raster()+
  geom_point(data = filter(Net_inds,patch==Fpatch,time==7000),aes(y=patch,x=time,fill=LD),size=5,pch=22, color="grey", stroke=2)+
  geom_point(data=contrast_data,aes(y=patch,x=time,fill=LD),size=5,pch=22, color="grey", stroke=2)+
  facet_grid(Community~Disp_text)+
  theme_bw(base_size = 12)+
  removeGrid()+
  scale_color_viridis(option = "D",name="Link\ndensity")+
  scale_fill_viridis(option = "D",name="Link\ndensity")+
  ylab("Patch")+
  xlab("Time")
ggsave(filename = "./Figures/Figure 4.png",width = 12,height = 8,dpi = 300)

#Figure S3####
Fig_S1a<-ggplot(filter(Net_turn_means_t,Part=="Gain",Scale=="Local network"),aes(x=Dispersal,y=WN_mean,color=Community, fill=Community))+
  #facet_grid(Part~Scale,scales="free")+
  geom_ribbon(aes(ymin=WN_lower,ymax=WN_upper),alpha=0.2, color=NA)+
  geom_line(size=1.2)+
  geom_line(data=filter(Net_turn_means3,Part=="Gain",Scale=="Local network"),aes(x=Dispersal,y=WN_mean),size=1,linetype=3)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_manual(values = ColV,labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
  scale_fill_manual(values = ColV,labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
  theme_bw(base_size = 12)+
  removeGrid()+
  ylab("Network dissimilarity (gains)")+
  xlab("Dispersal")+
  theme(legend.position="none")

Fig_S1b<-ggplot(filter(Net_turn_means_t,Part=="Loss",Scale=="Local network"),aes(x=Dispersal,y=WN_mean,color=Community, fill=Community))+
  #facet_grid(Part~Scale,scales="free")+
  geom_ribbon(aes(ymin=WN_lower,ymax=WN_upper),alpha=0.2, color=NA)+
  geom_line(size=1.2)+
  geom_line(data=filter(Net_turn_means3,Part=="Loss",Scale=="Local network"),aes(x=Dispersal,y=WN_mean),size=1,linetype=3)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_manual(values = ColV,labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
  scale_fill_manual(values = ColV,labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
  theme_bw(base_size = 12)+
  removeGrid()+
  ylab("Network dissimilarity (losses)")+
  xlab("Dispersal")+
  theme(legend.position=c(1,1),legend.justification=c(1,1))

Fig_S1c<-ggplot(filter(Trophic_shift_means,Part=="Gain",Scale=="Local network"),aes(x=Dispersal,y=WN_mean,color=Trophic,fill=Trophic))+
  #facet_grid(Part~Scale, scales="free")+
  geom_ribbon(aes(ymin=WN_lower,ymax=WN_upper),alpha=0.2, color=NA)+
  geom_line(size=1.2)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  theme_bw(base_size = 12)+
  removeGrid()+
  ylab("Network dissimilarity (gains)")+
  xlab("Dispersal")+
  theme(legend.position="none")

Fig_S1d<-ggplot(filter(Trophic_shift_means,Part=="Loss",Scale=="Local network"),aes(x=Dispersal,y=WN_mean,color=Trophic,fill=Trophic))+
  geom_ribbon(aes(ymin=WN_lower,ymax=WN_upper),alpha=0.2, color=NA)+
  geom_line(size=1.2)+
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

#Figures S4-S7####
ggplot(Net_inds,aes(y=patch,x=time,fill=N))+
  geom_raster()+
  geom_point(data = filter(Net_inds,patch==Fpatch,time==7000),aes(y=patch,x=time,fill=N),size=5,pch=22, color="grey", stroke=2)+
  geom_point(data=contrast_data,aes(y=patch,x=time,fill=N),size=5,pch=22, color="grey", stroke=2)+
  facet_grid(Community~Disp_text)+
  theme_bw(base_size = 12)+
  removeGrid()+
  scale_color_viridis(option = "D",name="# of species")+
  scale_fill_viridis(option = "D",name="# of species")+
  ylab("Patch")+
  xlab("Time")
ggsave(filename = "./Figures/Figure S4.png",width = 12,height = 8,dpi = 300)

ggplot(Net_inds,aes(y=patch,x=time,fill=Ltot))+
  geom_raster()+
  geom_point(data = filter(Net_inds,patch==Fpatch,time==7000),aes(y=patch,x=time,fill=Ltot),size=5,pch=22, color="grey", stroke=2)+
  geom_point(data=contrast_data,aes(y=patch,x=time,fill=Ltot),size=5,pch=22, color="grey", stroke=2)+
  facet_grid(Community~Disp_text)+
  theme_bw(base_size = 12)+
  removeGrid()+
  scale_color_viridis(option = "D",name="# of links")+
  scale_fill_viridis(option = "D",name="# of links")+
  ylab("Patch")+
  xlab("Time")
ggsave(filename = "./Figures/Figure S5.png",width = 12,height = 8,dpi = 300)

ggplot(Net_inds,aes(y=patch,x=time,fill=Nestedness))+
  geom_raster()+
  geom_point(data = filter(Net_inds,patch==Fpatch,time==7000),aes(y=patch,x=time,fill=Nestedness),size=5,pch=22, color="grey", stroke=2)+
  geom_point(data=contrast_data,aes(y=patch,x=time,fill=Nestedness),size=5,pch=22, color="grey", stroke=2)+
  facet_grid(Community~Disp_text)+
  theme_bw(base_size = 12)+
  removeGrid()+
  scale_color_viridis(option = "D",name="Nestedness")+
  scale_fill_viridis(option = "D",name="Nestedness")+
  ylab("Patch")+
  xlab("Time")
ggsave(filename = "./Figures/Figure S6.png",width = 12,height = 8,dpi = 300)

ggplot(Net_inds,aes(y=patch,x=time,fill=Cbar))+
  geom_raster()+
  geom_point(data = filter(Net_inds,patch==Fpatch,time==7000),aes(y=patch,x=time,fill=Cbar),size=5,pch=22, color="grey", stroke=2)+
  geom_point(data=contrast_data,aes(y=patch,x=time,fill=Cbar),size=5,pch=22, color="grey", stroke=2)+
  facet_grid(Community~Disp_text)+
  theme_bw(base_size = 12)+
  removeGrid()+
  scale_color_viridis(option = "D",name="Compartmentalization")+
  scale_fill_viridis(option = "D",name="Compartmentalization")+
  ylab("Patch")+
  xlab("Time")
ggsave(filename = "./Figures/Figure S7.png",width = 12,height = 8,dpi = 300)


#Fig S8####
Net_ind.df<-Net_ind.df%>%
  mutate(Trophic_levels=ifelse(Community!="Food web",NA,Trophic_levels))

Net_ind_long<-gather(Net_ind.df,key = Measure, value = Proportion,N:Trophic_levels)

Net_ind_long_means<-Net_ind_long%>%
  group_by(Dispersal,Community, Scale,Measure)%>%
  summarise_each(funs(mean(.,na.rm=T),lower=quantile(.,probs=c(0.25),na.rm=T),upper=quantile(.,probs=c(0.75),na.rm=T)))

Net_ind_long_means<-filter(Net_ind_long_means,Scale=="Local patch")
Net_ind_long_means<-filter(Net_ind_long_means,Measure=="N" | Measure=="Ltot" | Measure=="LD"| Measure=="Trophic_levels"|Measure=="Cbar"| Measure=="Nestedness")

Net_ind_long_means<-Net_ind_long_means%>%
  mutate(Measure=replace(Measure,Measure=="N","# of species"),
         Measure=replace(Measure,Measure=="Ltot","# of links"),
         Measure=replace(Measure,Measure=="LD","Link density"),
         Measure=replace(Measure,Measure=="Trophic_levels","Trophic levels"),
         Measure=replace(Measure,Measure=="Cbar","Compartmentalization"))


Net_ind_long_means$Measure<-factor(Net_ind_long_means$Measure,levels=c("# of species","# of links","Link density","Nestedness","Compartmentalization","Trophic levels"),ordered = T)


Net_ind_long_means_3<-Net_ind_long_means%>%
  filter(Community=="Food web",Dispersal>=0.05)

Net_ind_long_means<-Net_ind_long_means%>%
  mutate(Proportion_mean=replace(Proportion_mean,Community=="Food web" & Dispersal>0.05,NA))

ggplot(Net_ind_long_means,aes(x=Dispersal,y=Proportion_mean,color=Community, fill=Community))+
  facet_wrap(~Measure,scales="free_y")+
  geom_ribbon(aes(ymin=Proportion_lower,ymax=Proportion_upper),alpha=0.2, color=NA)+
  geom_line(size=1.2)+
  geom_line(data=Net_ind_long_means_3,aes(x=Dispersal,y=Proportion_mean),size=1,linetype=3)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c(0.0001,0.001,0.01,0.1,1))+
  scale_color_manual(values = ColV,labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
  scale_fill_manual(values = ColV,labels=c("No interactions", "Competition", "Mixed","Food web"),name="")+
  theme_bw(base_size = 12)+
  removeGrid()+
  geom_hline(yintercept = 1, linetype=2)+
  scale_y_continuous(breaks = c(0,0.5,1,1.5,2,2.5))+
  ylab("Propotion of initial network")+
  xlab("Dispersal")+
  theme(legend.position=c(1,1),legend.justification=c(1,1.6))
ggsave("./Figures/Figure S8.pdf",width = 10,height = 6,scale = 1)

#Fig S9####

int_hist.df<-rbind(int_hist_f(Com=Com_list$`Comp 1`,Ints = Int_list$BI,Fpatch = Fpatch,trophic = F,dispersal = 0.001,community = "Competition"),
      int_hist_f(Com=Com_list$`Comp 2`,Ints = Int_list$BI,Fpatch = Fpatch,trophic = F,dispersal = 0.01,community = "Competition"),
      int_hist_f(Com=Com_list$`Comp 3`,Ints = Int_list$BI,Fpatch = Fpatch,trophic = F,dispersal = 0.1,community = "Competition"),
      int_hist_f(Com=Com_list$`Comp 4`,Ints = Int_list$BI,Fpatch = Fpatch,trophic = F,dispersal = 0.5,community = "Competition"),
      int_hist_f(Com=Com_list$`Mixed 1`,Ints = Int_list$BM,Fpatch = Fpatch,trophic = F,dispersal = 0.001,community = "Mixed"),
      int_hist_f(Com=Com_list$`Mixed 2`,Ints = Int_list$BM,Fpatch = Fpatch,trophic = F,dispersal = 0.01,community = "Mixed"),
      int_hist_f(Com=Com_list$`Mixed 3`,Ints = Int_list$BM,Fpatch = Fpatch,trophic = F,dispersal = 0.1,community = "Mixed"),
      int_hist_f(Com=Com_list$`Mixed 4`,Ints = Int_list$BM,Fpatch = Fpatch,trophic = F,dispersal = 0.5,community = "Mixed"),
      int_hist_f(Com=Com_list$`FW 1`,Ints = Int_list$B3,Fpatch = Fpatch,trophic = T,dispersal = 0.001,community = "Food web"),
      int_hist_f(Com=Com_list$`FW 2`,Ints = Int_list$B3,Fpatch = Fpatch,trophic = T,dispersal = 0.01,community = "Food web"),
      int_hist_f(Com=Com_list$`FW 3`,Ints = Int_list$B3,Fpatch = Fpatch,trophic = T,dispersal = 0.1,community = "Food web"),
      int_hist_f(Com=Com_list$`FW 4`,Ints = Int_list$B3,Fpatch = Fpatch,trophic = T,dispersal = 0.5,community = "Food web"))

int_hist.df$Time<-factor(int_hist.df$Time,levels=c("Pre-change","Post-change"),ordered = T)

int_quants<-int_hist.df%>%
  group_by(Time,Dispersal,Community)%>%
  summarise(Quantile25=quantile(Real_ints,probs = c(0.25)),
            Quantile50=quantile(Real_ints,probs = c(0.5)),
            Quantile75=quantile(Real_ints,probs = c(0.75)))
int_quants<-gather(int_quants,key = Quantile,value = Value,Quantile25:Quantile75)
      

ggplot(int_hist.df,aes(x=Real_ints,color=Time, fill=Time))+
  geom_histogram(position="identity", bins=20)+
  theme(panel.margin=unit(3,"mm"))+
  facet_grid(Community~Dispersal, scales="free_y")+
  scale_fill_manual(values = c("grey",NA))+
  scale_color_brewer(palette = "Set1")+
  xlab("Realized interaction strength")+
  geom_vline(data=filter(int_quants,Quantile=="Quantile75"),aes(xintercept = Value, color=Time), linetype=2)
ggsave("./Figures/Figure S9.pdf",width = 12,height = 8)

int_scatter<-spread(data=int_hist.df,key = Time,value =  Real_ints)
int_scatter[is.na(int_scatter)]<-0
names(int_scatter)[4:5]<-c("pre","post")

ggplot(int_scatter,aes(x=pre,y=post))+
  geom_point()+
  facet_grid(Community~Dispersal)+
  geom_abline(intercept=0,slope=1)
ggsave("./Figures/Int scatter.pdf",width = 12,height = 8)

