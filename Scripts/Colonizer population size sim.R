#simulation to test the density dependence of colonization#
colonization_sim<-function(){
  nSpecies<-30
  
  weight=1/80*3
  
  dispV<-exp(seq(log(0.001),log(100),by=0.5))
  #dispV<-exp(seq(log(0.001),log(0.01),by=0.5))
  envV<-seq(2,9,by = 0.5)
  
  envStart<-10
  T_Opt<-c(1,rep(envStart,nSpecies-1))
  T_Norm<-apply(t(T_Opt),2,dnorm,sd=50,x=envStart)*300
  A<-(T_Norm-max(T_Norm))
  
  C<-rep(0.05,nSpecies)
  
  #competitive
  repeat{
    b11=-.15
    bdiag1=-.2
    BB=b11*matrix(runif(nSpecies*nSpecies),nSpecies,nSpecies)
    diag(BB)<-bdiag1
    BB=weight*BB
    BI<-BB
    BN=diag(diag(BB))
    
    BN1<-BN
    diag(BN1)<-0
    BI1<-BI
    diag(BI1)<-0
    
    if(solve(-BI,C+A)[1]<=0) {break}}
  
  
  #mixed
  repeat{
    BB<-matrix(-1,nSpecies,nSpecies)
    int.n<-sum(BB[upper.tri(BB)])*-1
    BB[upper.tri(BB)][sample(int.n, replace=F)<=(0.35*int.n)]<-0.5
    BB[lower.tri(BB)][t(BB)[lower.tri(BB)]>0][sample(0.35*int.n, replace=F)<(0.10*int.n)]<-0.5
    BM<-BB*-BI
    
    
    if(solve(-BM,C+A)[1]<=0) {break}}
  
  #Foodweb
  nprey<-5
  npred1<-7
  npred2<-11
  nSpeciesWeb<-sum(nprey,npred1,npred2)
  
  preyV<-1:nprey
  pred1V<-(nprey+1):(nprey+npred1)
  pred2V<-(nSpeciesWeb-npred2+1):(nSpeciesWeb)
  trophicV<-factor(c(rep("plant",nprey),rep("herbivore",npred1),rep("predator",npred2)),levels=c("plant","herbivore","predator"),ordered = T)
  
  C3<-c(rep(0.05,nprey),rep(0,nSpeciesWeb-nprey))
  
  T_Opt3<-rep(envStart,nSpeciesWeb)
  T_Opt3[match(unique(trophicV),trophicV)]<-1
  T_Norm<-apply(t(T_Opt3),2,dnorm,sd=50,x=envStart)*300
  A3<-(T_Norm-max(T_Norm))
  
  
  b11=-0.1
  b12=-0.3
  b21=0.1
  b23=-.1
  b32=.08
  bdiag1=-.2
  bdiag2=-.15
  
  
  #tritrophic BB Matrix####
  repeat{
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
    diag(BB[(nprey+npred1+1):nSpeciesWeb,(nprey+npred1+1):nSpeciesWeb])<-bdiag2
    BB=weight*BB
    
    B3<-BB
    if(solve(-B3,C3+A3)[1]<=0) {break}}
  
  
  for(d in 1:length(dispV)){
    print(d)
    for(env in 1:length(envV)){
      Env<-c(envStart,envV[env])
      
      XM<-X<-XI<-matrix(NA,nSpecies,2000)
      XM[,1]<-X[,1]<-XI[,1]<-10
      
      X3<-matrix(NA,nSpeciesWeb,2000)
      X3[,1]<-10
      
      
      T_Norm<-apply(t(T_Opt),2,dnorm,sd=50,x=Env)*300
      A<-(T_Norm-max(T_Norm))
      
      T_Norm<-apply(t(T_Opt3),2,dnorm,sd=50,x=Env)*300
      A3<-(T_Norm-max(T_Norm))
      
      for(l in 1:999){
        # X[,l+1]<-X[,l]*exp(C+BN%*%X[,l]+A[1,])
        # X[X<10^-3]<-0
        
        XI[,l+1]<-XI[,l]*exp(C+BI%*%XI[,l]+A[1,])
        XI[XI<10^-3]<-0
        
        XM[,l+1]<-XM[,l]*exp(C+BM%*%XM[,l]+A[1,])
        XM[XM<10^-3]<-0
        
        X3[,l+1]<-X3[,l]*exp(C3+B3%*%X3[,l]+A3[1,])
        X3[X3<10^-3]<-0
      }
      
      X[1,l+1]<-X[1,l+1]+dispV[d]
      XI[1,l+1]<-XI[1,l+1]+dispV[d]
      XM[1,l+1]<-XM[1,l+1]+dispV[d]
      X3[match(unique(trophicV),trophicV),l+1]<-X3[match(unique(trophicV),trophicV),l+1]+dispV[d]
      
      for(l in 1000:1999){
        # X[,l+1]<-X[,l]*exp(C+BN%*%X[,l]+A[2,])
        # X[X<10^-3]<-0
        
        XI[,l+1]<-XI[,l]*exp(C+BI%*%XI[,l]+A[2,])
        XI[XI<10^-3]<-0
        
        XM[,l+1]<-XM[,l]*exp(C+BM%*%XM[,l]+A[2,])
        XM[XM<10^-3]<-0
        
        X3[,l+1]<-X3[,l]*exp(C3+B3%*%X3[,l]+A3[2,])
        X3[X3<10^-3]<-0
      }
      hold<-data.frame(Env=envV[env],Dispersal=dispV[d],Community=c("Competitive","Mixed","Plant","Herbivore","Predator"),Rep=r,Invade=c(XI[1,l]>0,XM[1,l]>0,X3[match(unique(trophicV),trophicV),l]>0))
      if(env == 1 & d == 1){
        Results<-hold
      } else {
        Results<-rbind(Results,hold)
      }
    }
  }
  return(Results)
}

library(ggplot2)
library(dplyr)
library(ggExtra)
library(viridis)
library(RColorBrewer)
library(doParallel)
library(foreach)

cl<-makeCluster(2)
registerDoParallel(cl)
getDoParWorkers()

reps<-2

system.time(Sim_data_parallel<-foreach(r = 1:reps) %dopar% colonization_sim())

stopCluster(cl)

Sim_data<-do.call(rbind,Sim_data_parallel)

Result_means<-Sim_data%>%
  group_by(Dispersal,Env,Community)%>%
  summarise(Invade_prop=mean(Invade))

Result_means$Community<-factor(Result_means$Community,levels = c("Competitive","Mixed","Plant","Herbivore","Predator"),ordered = T)

ggplot(Result_means,aes(x=Dispersal,y=Env,fill=Invade_prop,z=Invade_prop))+
  facet_wrap(~Community)+
  geom_tile(color=NA) + 
  #stat_contour(aes(colour = ..level..))+
  scale_x_log10()+
  theme_bw()+
  xlab("Colonizing population size")+
  ylab("Environmental missmatch")+
  scale_fill_continuous(name="Colonization\nsuccess")
ggsave(filename = "./Figures/Colonization success.png",width = 8,height = 6)
