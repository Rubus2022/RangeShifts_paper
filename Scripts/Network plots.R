library(vegan)
library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(dplyr)
library(betalink)
library(igraph)
library(RInSp)
library(viridis)
library(NetIndices)

source("./Functions/plotting_functions.R")


dispV<-c(0.001,0.01,0.1,0.5)
reps<-1
dd<-c(0.3,0.2,0.1)#kernel decay strength
FoodWeb<-c("NoInt","Comp","Mixed","Plants","Herb","Pred") 

nCom<-200

Tmax<-3000
burnL<-2000
burn<-rep(1,burnL)
StressV<-c(burn,seq(1,Tmax),rep(Tmax,burnL))
maxStress<-20
Stress<-seq(0,maxStress,(maxStress)/(Tmax-1))
maxEnv1<-80
ComStart<-rev(seq(1,(maxEnv1),(maxEnv1-1)/(nCom-1)))

Disp1<-dist(seq(1:(nCom*2)))
Disp2<-matrix(0,nCom*2,nCom*2)
Disp2[lower.tri(Disp2)]<-Disp1;Disp2[upper.tri(Disp2)]<-rev(Disp1)
Disp<-exp(-dd[1]*Disp2)-diag(nrow(Disp2))
Disp<-Disp[((nCom/2)+1):((nCom/2)+nCom),((nCom/2)+1):((nCom/2)+nCom)]+rbind(apply(Disp[1:(nCom/2),((nCom/2)+1):((nCom/2)+nCom)],2,rev),apply(Disp[((nCom/2)+nCom+1):(nCom*2),((nCom/2)+1):((nCom/2)+nCom)],2,rev))

Disp_pl<-decostand(Disp,"total",2)

Disp<-exp(-dd[2]*Disp2)-diag(nrow(Disp2))
Disp<-Disp[((nCom/2)+1):((nCom/2)+nCom),((nCom/2)+1):((nCom/2)+nCom)]+rbind(apply(Disp[1:(nCom/2),((nCom/2)+1):((nCom/2)+nCom)],2,rev),apply(Disp[((nCom/2)+nCom+1):(nCom*2),((nCom/2)+1):((nCom/2)+nCom)],2,rev))
Disp_h<-decostand(Disp,"total",2)

Disp<-exp(-dd[3]*Disp2)-diag(nrow(Disp2))
Disp<-Disp[((nCom/2)+1):((nCom/2)+nCom),((nCom/2)+1):((nCom/2)+nCom)]+rbind(apply(Disp[1:(nCom/2),((nCom/2)+1):((nCom/2)+nCom)],2,rev),apply(Disp[((nCom/2)+nCom+1):(nCom*2),((nCom/2)+1):((nCom/2)+nCom)],2,rev))
Disp_pr<-decostand(Disp,"total",2)


#competitive
nprey<-80
npred1<-0
npred2<-0

n=nprey+npred1+npred2
weight=1/80*3

b11=-.15
bdiag1=-.2
BB=b11*matrix(runif(nprey*nprey),nprey,nprey)
diag(BB)<-bdiag1
BB=weight*BB
BI<-BB
BN=diag(diag(BB))

C<-rep(0.05,n)

T_Opt<-seq(1,maxEnv1,(maxEnv1-1)/(nprey-1))
T_Norm<-apply(t(T_Opt),2,dnorm,sd=50,x=seq(1,maxStress+maxEnv1))*300
A<-(T_Norm-max(T_Norm))

#mixed
BB<-matrix(-1,n,n)
int.n<-sum(BB[upper.tri(BB)])*-1
BB[upper.tri(BB)][sample(int.n, replace=F)<=(0.35*int.n)]<-0.5
BB[lower.tri(BB)][t(BB)[lower.tri(BB)]>0][sample(0.35*int.n, replace=F)<(0.10*int.n)]<-0.5
BM<-BB*-BI


#tri trophic
nprey<-40
npred1<-24
npred2<-16

preyV<-1:nprey
pred1V<-(nprey+1):(nprey+npred1)
pred2V<-(n-npred2+1):(n)

b11=-0.1
b12=-0.3
b21=0.1
b23=-.1
b32=.08
bdiag1=-.2
bdiag2=-.15


#tritrophic BB Matrix####
B11=b11*matrix(runif(nprey*nprey),nprey,nprey)
B12=b12*matrix(runif(nprey*npred1),nprey,npred1)
B13=matrix(0,nprey,npred2)
B21=b21*matrix(runif(npred1*nprey),npred1,nprey)
B22=matrix(0,npred1,npred1)
B23=b23*matrix(runif(npred1*npred2),npred1,npred2)
B31=matrix(0,npred2,nprey)
B32=b32*matrix(runif(npred2*npred1),npred2,npred1)
B33=matrix(0,npred2,npred2)
BB=rbind(cbind(B11 ,B12, B13),cbind(B21,B22, B23),cbind(B31, B32, B33))
diag(BB)<-bdiag1
diag(BB[(nprey+npred1+1):n,(nprey+npred1+1):n])<-bdiag2
BB=weight*BB

B3<-BB

BN1<-BN
diag(BN1)<-0
BI1<-BI
diag(BI1)<-0
BM1<-BM
diag(BM1)<-0
B31<-B3
diag(B31)<-0

C3<-c(rep(0.05,nprey),rep(0,n-nprey))

T_Opt3<-c(seq(1,maxEnv1,(maxEnv1-1)/(nprey-1)),seq(1,maxEnv1,(maxEnv1-1)/(npred1-1)),seq(1,maxEnv1,(maxEnv1-1)/(npred2-1)))
T_Norm<-apply(t(T_Opt3),2,dnorm,sd=50,x=seq(1,maxStress+maxEnv1))*300
A3<-(T_Norm-max(T_Norm))

Com_list<-list()
for(d in 1:length(dispV)){
  disp<-dispV[d]
  print(d)
  
  XI<-matrix(10,n,nCom)
  XM<-X3<-X<-XI
  XIt<-XMt<-X3t<-Xt<-XI
  
  sampleV<-seq(2000,7000,by=50)
  X_save<-array(NA,dim=c(n,nCom,length(sampleV)))
  XI_save<-XM_save<-X3_save<-X_save
  
  for(l in 1:(length(StressV)-1)){
    Xt<-X*exp(rep(C,nCom)+BN%*%X+t(A[(ComStart+Stress[StressV[l]]),]))+t(Disp_pl%*%t(X))*disp-disp*X
    Xt[Xt<10^-3]<-0
    X<-Xt
    
    XIt<-XI*exp(rep(C,nCom)+BI%*%XI+t(A[(ComStart+Stress[StressV[l]]),]))+t(Disp_pl%*%t(XI))*disp-disp*XI
    XIt[XIt<10^-3]<-0
    XI<-XIt
    
    XMt<-XM*exp(rep(C,nCom)+BM%*%XM+t(A[(ComStart+Stress[StressV[l]]),]))+t(Disp_pl%*%t(XM))*disp-disp*XM
    XMt[!is.finite(XMt)]<-0
    XMt[(XMt<10^-3)]<-0
    XM<-XMt
    
    X3t[preyV,]<-X3[preyV,]*exp(rep(C3[preyV],nCom)+B3[preyV,]%*%X3+t(A3[(ComStart+Stress[StressV[l]]),preyV]))+t(Disp_pl%*%t(X3[preyV,]))*disp-disp*X3[preyV,]
    X3t[pred1V,]<-X3[pred1V,]*exp(rep(C3[pred1V],nCom)+B3[pred1V,]%*%X3+t(A3[(ComStart+Stress[StressV[l]]),pred1V]))+t(Disp_h%*%t(X3[pred1V,]))*disp-disp*X3[pred1V,]
    X3t[pred2V,]<-X3[pred2V,]*exp(rep(C3[pred2V],nCom)+B3[pred2V,]%*%X3+t(A3[(ComStart+Stress[StressV[l]]),pred2V]))+t(Disp_pr%*%t(X3[pred2V,]))*disp-disp*X3[pred2V,]
    X3t[(X3t<10^-3)]<-0
    X3<-X3t
    if(l==2000){
      Xhold<-X
      XIhold<-XI
      XMhold<-XM
      X3hold<-X3
    }
    if(sum(l==(sampleV-1))==1){
      samp<-which(l==(sampleV-1))
      X_save[,,samp]<-X
      XI_save[,,samp]<-XI
      XM_save[,,samp]<-XM
      X3_save[,,samp]<-X3
    }
  }
  name<-paste("Comp",d)
  Com_list[[name]]<-XI_save
  name<-paste("Mixed",d)
  Com_list[[name]]<-XM_save
  name<-paste("FW",d)
  Com_list[[name]]<-X3_save
}

for(i in 1:length(dispV)){
  name<-paste("Mixed", i)
  Net_inds_1<-Net_ind_func(Com = Com_list[[name]],Ints = BM)
  Net_inds_1$Community<-"Mixed interactions"
  name<-paste("Comp", i)
  Net_inds_2<-Net_ind_func(Com = Com_list[[name]],Ints = BI)
  Net_inds_2$Community<-"Competition"
  name<-paste("FW", i)
  Net_inds_3<-Net_ind_func(Com = Com_list[[name]],Ints = B3)
  Net_inds_3$Community<-"Food web"
  hold<-rbind(Net_inds_1,Net_inds_2,Net_inds_3)
  hold$Dispersal<-dispV[i]
  if(i==1){
    Net_inds<-hold
  } else{ Net_inds<-rbind(Net_inds,hold)
  }
}

Net_inds$Community<-factor(Net_inds$Community,levels = c("Competition","Mixed interactions","Food web"),ordered = T)

sampleV<-seq(2000,7000,by=50)

Int_list<-list(BI=BI,BM=BM,B3=B3)

save(Net_inds,Com_list,Int_list,file="./Workspace/Range_shift_heatplots.RData")

ggplot(filter(Net_inds_3,Community!="No interactions"),aes(y=patch,x=time,fill=N))+
  geom_raster()+
  facet_grid(Community~Disp_text)+
  theme_bw(base_size = 12)+
  removeGrid()+
  #scale_colour_gradientn(colors=rainbow(100))+
  #scale_fill_gradientn(colors=rainbow(100,v = 1))+
  scale_color_viridis(option = "D")+
  scale_fill_viridis(option = "D")+
  ylab("Patch")+
  xlab("Time")
ggsave(filename = "./Figures/Species richness.png",width = 8,height = 8,dpi = 300)

ggplot(filter(Net_inds_3,Community!="No interactions"),aes(y=patch,x=time,fill=LD))+
  geom_raster()+
  facet_grid(Community~Disp_text)+
  theme_bw(base_size = 12)+
  removeGrid()+
  scale_color_viridis(option = "D")+
  scale_fill_viridis(option = "D")+
  ylab("Patch")+
  xlab("Time")
ggsave(filename = "./Figures/Link density.png",width = 8,height = 8,dpi = 300)

ggplot(filter(Net_inds_3,Community!="No interactions"),aes(y=patch,x=time,fill=C),color=NA)+
  geom_raster()+
  facet_grid(Community~Disp_text)+
  theme_bw(base_size = 12)+
  removeGrid()+
  scale_color_viridis(option = "D")+
  scale_fill_viridis(option = "D")+
  ylab("Patch")+
  xlab("Time")
ggsave(filename = "./Figures/Connectance.png",width = 8,height = 8,dpi = 300)

ggplot(filter(Net_inds_3,Community!="No interactions"),aes(y=patch,x=time,fill=Cbar),color=NA)+
  geom_raster()+
  facet_grid(Community~Disp_text)+
  theme_bw(base_size = 12)+
  removeGrid()+
  scale_color_viridis(option = "D")+
  scale_fill_viridis(option = "D")+
  ylab("Patch")+
  xlab("Time")
ggsave(filename = "./Figures/Compartmentalization.png",width = 8,height = 8,dpi = 300)

ggplot(filter(Net_inds_3,Community=="Food web"),aes(y=patch,x=time,fill=Trophic_levels),color=NA)+
  geom_raster()+
  facet_grid(.~Disp_text)+
  theme_bw(base_size = 12)+
  removeGrid()+
  scale_color_viridis(option = "D")+
  scale_fill_viridis(option = "D")+
  ylab("Patch")+
  xlab("Time")
ggsave(filename = "./Figures/Trophic levels.png",width = 8,height = 5,dpi = 300)

ggplot(filter(Net_inds_3,Community!="No interactions"),aes(y=patch,x=time,fill=Nestedness))+
  geom_raster()+
  facet_grid(Community~Disp_text)+
  theme_bw(base_size = 12)+
  removeGrid()+
  #scale_colour_gradientn(colors=rainbow(100))+
  #scale_fill_gradientn(colors=rainbow(100,v = 1))+
  scale_color_viridis(option = "D")+
  scale_fill_viridis(option = "D")+
  ylab("Patch")+
  xlab("Time")
ggsave(filename = "./Figures/Nestedness.png",width = 8,height = 8,dpi = 300)

