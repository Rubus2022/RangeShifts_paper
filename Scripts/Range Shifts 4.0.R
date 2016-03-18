library(vegan)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

source("./Functions/predict_err.r")

dispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1)
reps<-2
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

Bray.Curt<-array(NA,dim=c(reps,length(dispV),2,length(FoodWeb)),dimnames = list(1:reps,dispV,c("PA","AB"),FoodWeb))
Error_Rates<-array(NA,dim=c(reps,length(dispV),13,length(FoodWeb)),dimnames=list(seq(1:reps),dispV,c("Prevalence","Ov_Diag_Pow","Correct_class_rate","Sensitivity","Specificity","F_pos","F_neg","PPP", "NPP","Miss_class","Odds_ratio","Kappa","TSS"),FoodWeb))
Ext_Miss<-array(NA,dim=c(reps,length(dispV),length(FoodWeb)),dimnames = list(1:reps,dispV,FoodWeb))
Spatial.Insurance<-array(NA,dim=c(reps,length(dispV),8,length(FoodWeb)),dimnames = list(1:reps,dispV,c("R.SR_n","L.SR_n","L.CV","R.CV","R.SR","B.mass","B.mass.sd","L.SR"),FoodWeb))
Int_Stength<-array(NA,dim=c(80,reps,length(dispV),length(FoodWeb)),dimnames = list(1:80,1:reps,dispV,FoodWeb))
Cor_sd<-array(NA,dim=c(80,80,reps,length(dispV),length(FoodWeb)),dimnames=list(1:80,1:80,1:reps,dispV,FoodWeb))
Shift.Speed<-array(NA,dim=c(reps,length(dispV),length(FoodWeb)),dimnames=list(1:reps,dispV,FoodWeb))
Turnover<-array(NA,dim=c(reps,length(dispV),length(FoodWeb),5),dimnames=list(1:reps,dispV,FoodWeb,c("Pre","Post","Ext","Novel","Loc_ext")))

pb <- txtProgressBar(min = 0, max = reps, style = 3)
for(r in 1:reps){ 
  Sys.sleep(0.1)
  setTxtProgressBar(pb, r)
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
  
  
  for(d in 1:length(dispV)){
    disp<-dispV[d]
    print(d)
    
    XI=array(NA,dim=c(n,length(StressV),nCom))
    XI[,1,]=10
    XM<-X3<-X<-XI
    MeanInteract<-array(NA,dim=c(length(StressV)-1,n,4))
    
    for(l in 1:(length(StressV)-1)){
      X[,l+1,]<-X[,l,]*exp(rep(C,nCom)+BN%*%X[,l,]+t(A[(ComStart+Stress[StressV[l]]),]))+t(Disp_pl%*%t(X[,l,]))*disp-disp*X[,l,]
      X[,l+1,][(X[,l+1,]<10^-3)]<-0
      
      XI[,l+1,]<-XI[,l,]*exp(rep(C,nCom)+BI%*%XI[,l,]+t(A[(ComStart+Stress[StressV[l]]),]))+t(Disp_pl%*%t(XI[,l,]))*disp-disp*XI[,l,]
      XI[,l+1,][(XI[,l+1,]<10^-3)]<-0
      
      X3[preyV,l+1,]<-X3[preyV,l,]*exp(rep(C3[preyV],nCom)+B3[preyV,]%*%X3[,l,]+t(A3[(ComStart+Stress[StressV[l]]),preyV]))+t(Disp_pl%*%t(X3[preyV,l,]))*disp-disp*X3[preyV,l,]
      X3[preyV,l+1,][(X3[preyV,l+1,]<10^-3)]<-0
      X3[pred1V,l+1,]<-X3[pred1V,l,]*exp(rep(C3[pred1V],nCom)+B3[pred1V,]%*%X3[,l,]+t(A3[(ComStart+Stress[StressV[l]]),pred1V]))+t(Disp_h%*%t(X3[pred1V,l,]))*disp-disp*X3[pred1V,l,]
      X3[pred1V,l+1,][(X3[pred1V,l+1,]<10^-3)]<-0
      X3[pred2V,l+1,]<-X3[pred2V,l,]*exp(rep(C3[pred2V],nCom)+B3[pred2V,]%*%X3[,l,]+t(A3[(ComStart+Stress[StressV[l]]),pred2V]))+t(Disp_pr%*%t(X3[pred2V,l,]))*disp-disp*X3[pred2V,l,]
      X3[pred2V,l+1,][(X3[pred2V,l+1,]<10^-3)]<-0
      
      
      XM[,l+1,]<-XM[,l,]*exp(rep(C,nCom)+BM%*%XM[,l,]+t(A[(ComStart+Stress[StressV[l]]),]))+t(Disp_pl%*%t(XM[,l,]))*disp-disp*XM[,l,]
      XM[,l+1,][!is.finite(XM[,l+1,])]<-0
      XM[,l+1,][(XM[,l+1,]<10^-3)]<-0
      
      MeanInteract[l,,1]<-rowMeans(BN1%*%X[,l,])
      MeanInteract[l,,2]<-rowMeans(BI1%*%XI[,l,])
      MeanInteract[l,,3]<-rowMeans(BM1%*%XM[,l,])
      MeanInteract[l,,4]<-rowMeans(B31%*%X3[,l,])
      #print(l)
    }
    Turnover.df.temp<-rbind(Com_compare(X,Int_type = "None",Trophic = "No"),
                            Com_compare(XI,Int_type = "Competitive",Trophic = "No"),
                            Com_compare(XM,Int_type = "Mixed",Trophic = "No"),
                            Com_compare(X3[preyV,,],Int_type = "Plants",Trophic = "Yes"),
                            Com_compare(X3[pred1V,,],Int_type = "Herbivores",Trophic = "Yes"),
                            Com_compare(X3[pred2V,,],Int_type = "Predators",Trophic = "Yes"))
    
    
    Turnover.means.temp<-Turnover.df.temp%>%
      group_by(F_patch,I_patch,Interactions,Dispersal,Trophic)%>%
      mutate(Total_turnover=max(Turnover))%>%
      group_by(F_patch,Interactions,Dispersal,Trophic)%>%
      filter(Total_turnover==min(Total_turnover))%>%
      group_by(Type,Interactions,Dispersal,Trophic)%>%
      summarise(Turnover=mean(Turnover),Distance=mean(F_patch-I_patch))
    
    #calculate speed and variability of range shift
    rShift.df_temp<-rbind(rSpeedVary(X,Int_type = "None",Trophic = "No"),
                          rSpeedVary(XI,Int_type = "Competitive",Trophic = "No"),
                          rSpeedVary(XM,Int_type = "Mixed",Trophic = "No"),
                          rSpeedVary(X3[preyV,,],Int_type = "Plants",Trophic = "Yes"),
                          rSpeedVary(X3[pred1V,,],Int_type = "Herbivores",Trophic = "Yes"),
                          rSpeedVary(X3[pred2V,,],Int_type = "Predators",Trophic = "Yes"))
    
    
    if(r==1 & d==1){
      rShift.df<-rShift.df_temp
      Turn_mean<-Turnover.means.temp
    } else {rShift.df<-rbind(rShift.df,rShift.df_temp)
    Turn_mean<-rbind(Turn_mean,Turnover.means.temp)}
    
  }  
};close(pb)

save(dispV,Turn_mean,rShift.df,file = "Species Interactions.RData")

Turn_means<-Turn_mean%>%
  group_by(Type,Interactions,Dispersal,Trophic)%>%
  summarise(Turnover=mean(Turnover),Distance=mean(Distance))

ggplot(Turn_means,aes(x=Dispersal,y=Turnover,color=Interactions))+
  geom_line()+
  facet_grid(Type~Trophic,scale="free_y")+
  scale_x_log10()

ggplot(filter(Turn_means,Type=="Total"),aes(x=Dispersal,y=Distance,color=Interactions))+
  geom_line()+
  facet_grid(.~Trophic,scale="free_y")+
  scale_x_log10()


#Figures####
setwd("/Users/patrickthompson/Dropbox/Patrick/Species Interactions/Manuscript")
ColV<-c(8,1,"dodgerblue",brewer.pal(7,"Set1")[c(3,4,1)])
LtyV<-c(1,1,1,1,1,1)
ColT<-c(1,brewer.pal(3,"Dark2"))
lab.size<-1.4
axis.size<-1.2
pchV<-c(20,20,20,20,20,20)

####Paper Figures####
#Fig 2####
pdf("Figure 2.pdf",width=9,height=9)
par(mfrow=c(2,2),las=1, pty='s', oma=c(2,2,2,2),mar=c(3,5,1,2))
plot(dispV,colMeans(Bray.Curt[,,"AB","NoInt"]),lwd=3, log="x", type='n', ylim=c(0,1.05), ylab="Bray-Curtis similarity",col=ColV[1],xlab="Dispersal", xaxt='n',cex.lab=lab.size,cex.axis=axis.size)
for(i in 1:3){
  lines(dispV,colMeans(Bray.Curt[,,"AB",i],na.rm=T),lwd=3, col=ColV[i], lty=LtyV[i])
  plotCI(dispV,colMeans(Bray.Curt[,,"AB",i],na.rm=T),uiw=apply(Bray.Curt[,,"AB",i],2,sd,na.rm=T), lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
mtext("a",3,line=0.25,cex=1.5,adj=-.21)

plot(dispV,colMeans(Bray.Curt[,,"AB","NoInt"]),lwd=3, log="x", type='n', ylim=c(0,1.05), ylab="Bray-Curtis similarity",col=ColV[1],xlab="Dispersal", xaxt='n',cex.lab=lab.size,cex.axis=axis.size)
for(i in 4:6){
  lines(dispV,colMeans(Bray.Curt[,,"AB",i],na.rm=T),lwd=3, col=ColV[i], lty=LtyV[i])
  plotCI(dispV,colMeans(Bray.Curt[,,"AB",i],na.rm=T),uiw=apply(Bray.Curt[,,"AB",i],2,sd,na.rm=T), lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
mtext("b",3,line=0.25,cex=1.5,adj=-.21)

plot(colMeans(Shift.Speed[,,"Comp"],na.rm=T)~dispV,log='x', type='n', ylim=c(0,75), ylab=NA, xlab="Dispersal",xaxt='n',cex.lab=lab.size,cex.axis=axis.size)
title(ylab="Range shift rate variation\n(interspecific sd)",mgp=c(2.5,1,0),cex.lab=lab.size)
for(i in 1:3){
  lines(colMeans(Shift.Speed[,,i],na.rm=T)~dispV, lwd=3,col=ColV[i], lty=LtyV[i])
  plotCI(dispV,colMeans(Shift.Speed[,,i],na.rm=T),uiw=apply(Shift.Speed[,,i],2,sd,na.rm=T), lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
mtext("c",3,line=0.25,cex=1.5,adj=-.21)
legend("topright",c("No interactions","Competitive","Mixed"),lwd=3,lty=1,col=ColV[1:3],pch=pchV[1:3], bty='n', cex=1.2)

plot(colMeans(Shift.Speed[,,"Comp"],na.rm=T)~dispV,log='x', type='n', ylim=c(0,75), ylab=NA, xlab="Dispersal",xaxt='n',cex.lab=lab.size,cex.axis=axis.size)
title(ylab="Range shift rate variation\n(interspecific sd)",mgp=c(2.5,1,0),cex.lab=lab.size)
for(i in 4:6){
  lines(colMeans(Shift.Speed[,,i],na.rm=T)~dispV, lwd=3,col=ColV[i], lty=LtyV[i])
  plotCI(dispV,colMeans(Shift.Speed[,,i],na.rm=T),uiw=apply(Shift.Speed[,,i],2,sd,na.rm=T), lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
mtext("d",3,line=0.25,cex=1.5,adj=-.21)
legend("topright",inset=c(0,0),c("Plants","Herbivores","Carnivores"),lwd=3,lty=1,col=ColV[4:6],pch=pchV[4:6], bty='n', cex=1.2)
dev.off()

pdf("Figure 2_PA.pdf",width=9,height=9)
par(mfrow=c(2,2),las=1, pty='s', oma=c(2,2,2,2),mar=c(3,5,1,2))
plot(dispV,colMeans(Bray.Curt[,,"PA","NoInt"]),lwd=3, log="x", type='n', ylim=c(0,1.05), ylab="Bray-Curtis similarity",col=ColV[1],xlab="Dispersal", xaxt='n',cex.lab=lab.size,cex.axis=axis.size)
for(i in 1:3){
  lines(dispV,colMeans(Bray.Curt[,,"PA",i],na.rm=T),lwd=3, col=ColV[i], lty=LtyV[i])
  plotCI(dispV,colMeans(Bray.Curt[,,"PA",i],na.rm=T),uiw=apply(Bray.Curt[,,"PA",i],2,sd,na.rm=T), lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
mtext("a",3,line=0.25,cex=1.5,adj=-.21)

plot(dispV,colMeans(Bray.Curt[,,"PA","NoInt"]),lwd=3, log="x", type='n', ylim=c(0,1.05), ylab="Bray-Curtis similarity",col=ColV[1],xlab="Dispersal", xaxt='n',cex.lab=lab.size,cex.axis=axis.size)
for(i in 4:6){
  lines(dispV,colMeans(Bray.Curt[,,"PA",i],na.rm=T),lwd=3, col=ColV[i], lty=LtyV[i])
  plotCI(dispV,colMeans(Bray.Curt[,,"PA",i],na.rm=T),uiw=apply(Bray.Curt[,,"PA",i],2,sd,na.rm=T), lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
mtext("b",3,line=0.25,cex=1.5,adj=-.21)

plot(colMeans(Shift.Speed[,,"Comp"],na.rm=T)~dispV,log='x', type='n', ylim=c(0,75), ylab=NA, xlab="Dispersal",xaxt='n',cex.lab=lab.size,cex.axis=axis.size)
title(ylab="Range shift rate variation\n(interspecific sd)",mgp=c(2.5,1,0),cex.lab=lab.size)
for(i in 1:3){
  lines(colMeans(Shift.Speed[,,i],na.rm=T)~dispV, lwd=3,col=ColV[i], lty=LtyV[i])
  plotCI(dispV,colMeans(Shift.Speed[,,i],na.rm=T),uiw=apply(Shift.Speed[,,i],2,sd,na.rm=T), lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
mtext("c",3,line=0.25,cex=1.5,adj=-.21)
legend("topright",c("No interactions","Competitive","Mixed"),lwd=3,lty=1,col=ColV[1:3],pch=pchV[1:3], bty='n', cex=1.2)

plot(colMeans(Shift.Speed[,,"Comp"],na.rm=T)~dispV,log='x', type='n', ylim=c(0,75), ylab=NA, xlab="Dispersal",xaxt='n',cex.lab=lab.size,cex.axis=axis.size)
title(ylab="Range shift rate variation\n(interspecific sd)",mgp=c(2.5,1,0),cex.lab=lab.size)
for(i in 4:6){
  lines(colMeans(Shift.Speed[,,i],na.rm=T)~dispV, lwd=3,col=ColV[i], lty=LtyV[i])
  plotCI(dispV,colMeans(Shift.Speed[,,i],na.rm=T),uiw=apply(Shift.Speed[,,i],2,sd,na.rm=T), lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
mtext("d",3,line=0.25,cex=1.5,adj=-.21)
legend("topright",inset=c(0,0),c("Plants","Herbivores","Carnivores"),lwd=3,lty=1,col=ColV[4:6],pch=pchV[4:6], bty='n', cex=1.2)
dev.off()

#Fig 4####
pdf("Figure 4.pdf",width=9,height=9)
par(mfrow=c(2,2),las=1, pty='s', oma=c(2,2,2,2),mar=c(3,5,1,2))
plot(dispV,colMeans(Ext_Miss[,,2]),lwd=3, log="x", type='n', ylim=c(0,17),ylab="Extinction order mismatch",col=ColV[1],xlab="Dispersal", xaxt='n', cex.lab=lab.size,cex.axis=axis.size)
for(i in 1:3){
  lines(dispV,colMeans(Ext_Miss[,,i],na.rm=T),lwd=3, col=ColV[i],lty=LtyV[i])
  plotCI(dispV,colMeans(Ext_Miss[,,i],na.rm=T),uiw=apply(Ext_Miss[,,i],2,sd,na.rm=T), lwd=2, add=T, pch=pchV, col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
legend("topright",c("No interactions","Competitive","Mixed"),lwd=3,lty=1,col=ColV[1:3],pch=pchV[1:3], bty='n', cex=1.2)
mtext("a",3,line=0.25,cex=1.5,adj=-.21)

plot(dispV,colMeans(Ext_Miss[,,2]),lwd=3, log="x", type='n', ylim=c(0,17),ylab="Extinction order mismatch",col=ColV[1],xlab="Dispersal", xaxt='n', cex.lab=lab.size,cex.axis=axis.size)
for(i in 4:6){
  lines(dispV,colMeans(Ext_Miss[,,i],na.rm=T),lwd=3, col=ColV[i],lty=LtyV[i])
  plotCI(dispV,colMeans(Ext_Miss[,,i],na.rm=T),uiw=apply(Ext_Miss[,,i],2,sd,na.rm=T), lwd=2, add=T, pch=pchV, col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
legend("topright",c("Plants","Herbivores","Carnivores"),lwd=3,lty=1,col=ColV[4:6],pch=pchV[4:6], bty='n', cex=1.2)
mtext("b",3,line=0.25,cex=1.5,adj=-.21)
plot(NA)
plot(NA)
dev.off()

# intersp correlation sd####
pdf("Interspecific correlation sd.pdf",height=12,width=6)
par(mfrow=c(2,1))
plot(apply(Cor_sd[,,1,,"Mixed"],3,mean,na.rm=T)~dispV,log='x', type='n', ylim=c(0.02,0.1), ylab="Temporal variation in interspecific correlations", xlab="Dispersal")
for(i in 1:6){
  lines(dispV,apply(Cor_sd[,,1,,i],3,mean,na.rm=T),lwd=3, col=ColV[i])
}
legend("bottomleft", legend=c("Competition","Mixed","Trophic"), lwd=3,col=c(1,4,2), bty='n')

plot(apply(Cor_sd[,,1,,"Arb"],3,mean,na.rm=T)~dispV,log='x', type='n', ylim=c(0.02,0.06), ylab="Temporal variation in interspecific correlations", xlab="Dispersal")
lines(dispV,apply(Cor_sd[1:nprey,,1,,3],3,mean,na.rm=T),lwd=3, col=ColT[2])
lines(dispV,apply(Cor_sd[(nprey+1):(nprey+npred1),,1,,3],3,mean,na.rm=T),lwd=3, col=ColT[3])
lines(dispV,apply(Cor_sd[(n-npred2):(n),,1,,3],3,mean,na.rm=T),lwd=3, col=ColT[4])
legend("bottomleft", bty='n', col=ColT[2:4],c("Plants","Herbivores","Predators"),lwd=3)
dev.off()

#plot interaction strength####
pdf("Interaction strength.pdf",height=6,width=6)
par(mfrow=c(1,1),mar=c(5,5,5,5), oma=c(0,0,0,0))
plot(apply(Int_Stength[,,,"Comp"],3,mean)~dispV,log='x', type='l',col=ColV[2], ylim=c(-0.04,0.02),lwd=3, ylab="Biotic interaction strength", xlab="Dispersal",xaxt='n')
lines(apply(Int_Stength[preyV,,,"Plants"],3,mean)~dispV, col=ColV[4],lwd=3)
lines(apply(Int_Stength[pred1V,,,"Herb"],3,mean)~dispV, col=ColV[5],lwd=3)
lines(apply(Int_Stength[pred2V,,,"Pred"],3,mean)~dispV, col=ColV[6],lwd=3)
lines(apply(Int_Stength[,,,"Mixed"],3,mean)~dispV, col=ColV[3],lwd=3)
lines(apply(Int_Stength[,,,"NoInt"],3,mean)~dispV, col=ColV[1],lwd=3)
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
dev.off()

#plot variation in speed of range shift####
pdf("Range shift speed sd.pdf",height=6,width=6)
plot(colMeans(Shift.Speed[,,"Comp",1],na.rm=T)~dispV,log='x', type='n', ylim=c(0,50), ylab="Range shift rate (interspecific sd)", xlab="Dispersal")
for(i in 1:5){
  lines(colMeans(Shift.Speed[,,i,1],na.rm=T)~dispV, lwd=3,col=ColV[i])}
dev.off()

#Supplemental####
pdf("Figure S1.pdf", height =9, width=9)
par(mfrow=c(2,2),las=1, pty='s', oma=c(2,2,2,2),mar=c(3,5,1,2))
plot(dispV,colMeans(Error_Rates[,,"Correct_class_rate","NoInt"],na.rm=T),lwd=3, log="x", type='n', ylim=c(0,1), ylab="Correct classification rate",col=ColV[1],xlab="Dispersal", xaxt='n',cex.lab=lab.size,cex.axis=axis.size)
for(i in 1:3){
  lines(dispV,colMeans(Error_Rates[,,"Correct_class_rate",i],na.rm=T),lwd=3, col=ColV[i],lty=LtyV[i])
  plotCI(dispV,colMeans(Error_Rates[,,"Correct_class_rate",i],na.rm=T),uiw=apply(Error_Rates[,,"Correct_class_rate",i],2,sd,na.rm=T), lwd=2, add=T, col=ColV[i], pch=pchV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
mtext("a",3,line=0.25,cex=1.5,adj=-.21)

plot(dispV,colMeans(Error_Rates[,,"Correct_class_rate","NoInt"],na.rm=T),lwd=3, log="x", type='n', ylim=c(0,1), ylab="Correct classification rate",col=ColV[1],xlab="Dispersal", xaxt='n',cex.lab=lab.size,cex.axis=axis.size)
for(i in 4:6){
  lines(dispV,colMeans(Error_Rates[,,"Correct_class_rate",i],na.rm=T),lwd=3, col=ColV[i],lty=LtyV[i])
  plotCI(dispV,colMeans(Error_Rates[,,"Correct_class_rate",i],na.rm=T),uiw=apply(Error_Rates[,,"Correct_class_rate",i],2,sd,na.rm=T), lwd=2, add=T, col=ColV[i], pch=pchV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
mtext("b",3,line=0.25,cex=1.5,adj=-.21)

plot(dispV,colMeans(Error_Rates[,,"TSS","NoInt"],na.rm=T),lwd=3, log="x", type='n', ylim=c(-1,1), ylab="True skill statistic",col=ColV[1],xlab="Dispersal", xaxt='n',cex.lab=lab.size,cex.axis=axis.size)
for(i in 1:3){
  lines(dispV,colMeans(Error_Rates[,,"TSS",i],na.rm=T),lwd=3, col=ColV[i],lty=LtyV[i])
  plotCI(dispV,colMeans(Error_Rates[,,"TSS",i],na.rm=T),uiw=apply(Error_Rates[,,"TSS",i],2,sd,na.rm=T), lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
mtext("c",3,line=0.25,cex=1.5,adj=-.21)
legend("bottomright",c("No interactions","Competitive","Mixed"),lwd=3,lty=1,col=ColV[1:3],pch=pchV[1:3], bty='n', cex=1.2)

plot(dispV,colMeans(Error_Rates[,,"TSS","NoInt"],na.rm=T),lwd=3, log="x", type='n', ylim=c(-1,1), ylab="True skill statistic",col=ColV[1],xlab="Dispersal", xaxt='n',cex.lab=lab.size,cex.axis=axis.size)
for(i in 4:6){
  lines(dispV,colMeans(Error_Rates[,,"TSS",i],na.rm=T),lwd=3, col=ColV[i],lty=LtyV[i])
  plotCI(dispV,colMeans(Error_Rates[,,"TSS",i],na.rm=T),uiw=apply(Error_Rates[,,"TSS",i],2,sd,na.rm=T), lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
mtext("d",3,line=0.25,cex=1.5,adj=-.21)
legend("bottomright",c("Plants","Herbivores","Carnivores"),lwd=3,lty=1,col=ColV[4:6],pch=pchV[4:6], bty='n', cex=1.2)
dev.off()

#Spatial insurance####
pdf("Lanscape spatial insurance.pdf",width=12,height=7.5)
par(mfrow=c(2,3),las=1, pty='s', oma=c(2,2,2,2),mar=c(5,5,1,2))
plot(dispV,colMeans(Spatial.Insurance[,,"L.SR",2]),lwd=3, log="x", type='n',ylim=c(0,1),ylab="Propotion of initial local SR",col=ColV[1],xlab="Dispersal", xaxt='n', cex.lab=lab.size)
for(i in 1:3){
  lines(dispV,colMeans(Spatial.Insurance[,,"L.SR",i],na.rm=T),lwd=3, col=ColV[i])
  plotCI(dispV,colMeans(Spatial.Insurance[,,"L.SR",i],na.rm=T),uiw=apply(Spatial.Insurance[,,"L.SR",i],2,sd,na.rm=T), lwd=2, add=T, pch=20, col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))
legend("bottomright",c("No interactions","Competitive","Mixed"),lwd=3,lty=1,col=ColV[1:3],pch=pchV[1:3], bty='n', cex=1.2)
mtext("a",3,line=0.25,cex=1.5,adj=-.21)

plot(dispV,colMeans(Spatial.Insurance[,,"R.SR",2]),lwd=3, log="x", type='n',ylim=c(0,1),ylab="Propotion of initial regional SR",col=ColV[1],xlab="Dispersal", xaxt='n', cex.lab=lab.size)
for(i in 1:3){
  lines(dispV,colMeans(Spatial.Insurance[,,"R.SR",i],na.rm=T),lwd=3, col=ColV[i])
  plotCI(dispV,colMeans(Spatial.Insurance[,,"R.SR",i],na.rm=T),uiw=apply(Spatial.Insurance[,,"R.SR",i],2,sd,na.rm=T), lwd=2, add=T, pch=20, col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))
mtext("b",3,line=0.25,cex=1.5,adj=-.21)

plot(dispV,colMeans(Spatial.Insurance[,,"B.mass",2]),lwd=3, log="x", type='n',ylim=c(0,1.2),ylab="Propotion of initial biomass",col=ColV[1],xlab="Dispersal", xaxt='n', cex.lab=lab.size)
for(i in 1:3){
  lines(dispV,colMeans(Spatial.Insurance[,,"B.mass",i],na.rm=T),lwd=3, col=ColV[i])
  plotCI(dispV,colMeans(Spatial.Insurance[,,"B.mass",i],na.rm=T),uiw=apply(Spatial.Insurance[,,"B.mass",i],2,sd,na.rm=T), lwd=2, add=T, pch=20, col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))
mtext("c",3,line=0.25,cex=1.5,adj=-.21)

plot(dispV,colMeans(Spatial.Insurance[,,"L.CV",2]),lwd=3, log="x", type='n',ylim=c(0,6),ylab="Local biomass CV",col=ColV[1],xlab="Dispersal", xaxt='n', cex.lab=lab.size)
for(i in 1:3){
  lines(dispV,colMeans(Spatial.Insurance[,,"L.CV",i],na.rm=T),lwd=3, col=ColV[i])
  plotCI(dispV,colMeans(Spatial.Insurance[,,"L.CV",i],na.rm=T),uiw=apply(Spatial.Insurance[,,"L.CV",i],2,sd,na.rm=T), lwd=2, add=T, pch=20, col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))
mtext("d",3,line=0.25,cex=1.5,adj=-.21)

plot(dispV,colMeans(Spatial.Insurance[,,"R.CV",2]),lwd=3, log="x", type='n',ylim=c(0,0.8),ylab="Regional biomass CV",col=ColV[1],xlab="Dispersal", xaxt='n', cex.lab=lab.size)
for(i in 1:3){
  lines(dispV,colMeans(Spatial.Insurance[,,"R.CV",i],na.rm=T),lwd=3, col=ColV[i])
  plotCI(dispV,colMeans(Spatial.Insurance[,,"R.CV",i],na.rm=T),uiw=apply(Spatial.Insurance[,,"R.CV",i],2,sd,na.rm=T), lwd=2, add=T, pch=20, col=ColV[i])
}
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))
mtext("e",3,line=0.25,cex=1.5,adj=-.21)
dev.off()


#Turnover plot####
Turn_prop<-array(NA,dim=c(reps,length(dispV),length(FoodWeb),4),dimnames = list(1:reps,dispV,FoodWeb,c("Turnover","Novel","Local ext","Ext")))

for(i in 1:length(FoodWeb)){
  Turn_prop[,,i,"Turnover"]<-(Turnover[,,i,"Novel"]+Turnover[,,i,"Ext"])/(Turnover[,,i,"Pre"]+Turnover[,,i,"Post"])
  Turn_prop[,,i,"Ext"]<-(Turnover[,,i,"Ext"]-Turnover[,,i,"Loc_ext"])/(Turnover[,,i,"Pre"]+Turnover[,,i,"Post"])
  Turn_prop[,,i,"Local ext"]<-(Turnover[,,i,"Loc_ext"])/(Turnover[,,i,"Pre"]+Turnover[,,i,"Post"])
  Turn_prop[,,i,"Novel"]<-(Turnover[,,i,"Novel"])/(Turnover[,,i,"Pre"]+Turnover[,,i,"Post"])
}

pdf("Assemblage turnover.pdf",width=8.5,height=11)
par(mfrow=c(3,2),las=1, pty='s', oma=c(3,3,3,3),mar=c(5,5,1,2))
plot(dispV,colMeans(Turn_prop[,,1,"Turnover"],na.rm=T)*100,lwd=3, col=ColV[i], type='n', ylim=c(0,100),log='x',xaxt='n',ylab="Species assemblage\nturnover (%)",xlab="Dispersal",cex.lab=lab.size,cex.axis=axis.size)
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
for(i in 1:3){
  lines(dispV,colMeans(Turn_prop[,,i,"Turnover"],na.rm=T)*100,lwd=3, col=ColV[i], lty=1)
  plotCI(dispV,colMeans(Turn_prop[,,i,"Turnover"],na.rm=T)*100,uiw=apply(Turn_prop[,,i,"Turnover"],2,sd,na.rm=T)*100, lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
mtext("a",3,line=0.25,cex=1.5,adj=-.21)
legend("topright",c("No interactions","Competitive","Mixed"),lwd=3,lty=1,col=ColV[1:3],pch=pchV[1:3], bty='n', cex=1.2)

plot(dispV,colMeans(Turn_prop[,,1,"Turnover"],na.rm=T)*100,lwd=3, col=ColV[i], type='n', ylim=c(0,100),log='x',xaxt='n',ylab="",xlab="Dispersal",cex.lab=lab.size,cex.axis=axis.size)
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
for(i in 4:6){
  lines(dispV,colMeans(Turn_prop[,,i,"Turnover"],na.rm=T)*100,lwd=3, col=ColV[i], lty=1)
  plotCI(dispV,colMeans(Turn_prop[,,i,"Turnover"],na.rm=T)*100,uiw=apply(Turn_prop[,,i,"Turnover"],2,sd,na.rm=T)*100, lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
legend("bottomleft",c("Plants","Herbivores","Carnivores"),lwd=3,lty=1,col=ColV[4:6],pch=pchV[4:6], bty='n', cex=1.2)
mtext("b",3,line=0.25,cex=1.5,adj=-.21)

dev.off()


pdf("Turnover.pdf",width=8.5,height=11)
par(mfrow=c(3,2),las=1, pty='s', oma=c(3,3,3,3),mar=c(5,5,1,2))
plot(dispV,colMeans(Turn_prop[,,1,"Ext"],na.rm=T)*100,lwd=3, col=ColV[i], type='n', ylim=c(0,100),log='x',xaxt='n',ylab="Assemblage turnover (%)\nfrom extinction",xlab="Dispersal",cex.lab=lab.size,cex.axis=axis.size)
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
for(i in 1:3){
  lines(dispV,colMeans(Turn_prop[,,i,"Ext"],na.rm=T)*100,lwd=3, col=ColV[i], lty=1)
  plotCI(dispV,colMeans(Turn_prop[,,i,"Ext"],na.rm=T)*100,uiw=apply(Turn_prop[,,i,"Ext"],2,sd,na.rm=T)*100, lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
mtext("a",3,line=0.25,cex=1.5,adj=-.21)
plot(dispV,colMeans(Turn_prop[,,i,"Ext"],na.rm=T)*100,lwd=3, col=ColV[i], lty=LtyV[i], type='n', ylim=c(0,100),log='x',xaxt='n',ylab="",xlab="Dispersal",cex.lab=lab.size,cex.axis=axis.size)
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
for(i in 4:6){
  lines(dispV,colMeans(Turn_prop[,,i,"Ext"],na.rm=T)*100,lwd=3, col=ColV[i], lty=1)
  plotCI(dispV,colMeans(Turn_prop[,,i,"Ext"],na.rm=T)*100,uiw=apply(Turn_prop[,,i,"Ext"],2,sd,na.rm=T)*100, lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
mtext("b",3,line=0.25,cex=1.5,adj=-.21)
plot(dispV,colMeans(Turn_prop[,,1,"Novel"],na.rm=T)*100,lwd=3, col=ColV[i], type='n', ylim=c(0,11),log='x',xaxt='n',ylab="Assemblage turnover (%)\nfrom novel species",xlab="Dispersal",cex.lab=lab.size,cex.axis=axis.size)
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
for(i in 1:3){
  lines(dispV,colMeans(Turn_prop[,,i,"Novel"],na.rm=T)*100,lwd=3, col=ColV[i], lty=1)
  plotCI(dispV,colMeans(Turn_prop[,,i,"Novel"],na.rm=T)*100,uiw=apply(Turn_prop[,,i,"Novel"],2,sd,na.rm=T)*100, lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
mtext("c",3,line=0.25,cex=1.5,adj=-.21)
plot(dispV,colMeans(Turn_prop[,,i,"Novel"],na.rm=T)*100,lwd=3, col=ColV[i], lty=LtyV[i], type='n', ylim=c(0,11),log='x',xaxt='n',ylab="",xlab="Dispersal",cex.lab=lab.size,cex.axis=axis.size)
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
for(i in 4:6){
  lines(dispV,colMeans(Turn_prop[,,i,"Novel"],na.rm=T)*100,lwd=3, col=ColV[i], lty=1)
  plotCI(dispV,colMeans(Turn_prop[,,i,"Novel"],na.rm=T)*100,uiw=apply(Turn_prop[,,i,"Novel"],2,sd,na.rm=T)*100, lwd=2, add=T, pch=pchV[i], col=ColV[i])
} 
mtext("d",3,line=0.25,cex=1.5,adj=-.21)
plot(dispV,colMeans(Turn_prop[,,1,"Local ext"],na.rm=T)*100,lwd=3, col=ColV[i], type='n', ylim=c(0,36),log='x',xaxt='n',ylab="Assemblage turnover (%)\nfrom local extirpations",xlab="Dispersal",cex.lab=lab.size,cex.axis=axis.size)
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
for(i in 1:3){
  lines(dispV,colMeans(Turn_prop[,,i,"Local ext"],na.rm=T)*100,lwd=3, col=ColV[i], lty=1)
  plotCI(dispV,colMeans(Turn_prop[,,i,"Local ext"],na.rm=T)*100,uiw=apply(Turn_prop[,,i,"Local ext"],2,sd,na.rm=T)*100, lwd=2, add=T, pch=pchV[i], col=ColV[i])
}
legend("topright",c("No interactions","Competitive","Mixed"),lwd=3,lty=1,col=ColV[1:3],pch=pchV[1:3], bty='n', cex=1.2)
mtext("e",3,line=0.25,cex=1.5,adj=-.21)
plot(dispV,colMeans(Turn_prop[,,i,"Local ext"],na.rm=T)*100,lwd=3, col=ColV[i], lty=LtyV[i], type='n', ylim=c(0,36),log='x',xaxt='n',ylab="",xlab="Dispersal",cex.lab=lab.size,cex.axis=axis.size)
axis(1,at=c(seq(0.0002,0.0009,0.0001),seq(0.002,0.009,0.001),seq(0.02,0.09,0.01),seq(0.2,0.9,0.1)),tcl=-0.3,labels=F)
axis(1,at=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"),cex.axis=axis.size)
for(i in 4:6){
  lines(dispV,colMeans(Turn_prop[,,i,"Local ext"],na.rm=T)*100,lwd=3, col=ColV[i], lty=1)
  plotCI(dispV,colMeans(Turn_prop[,,i,"Local ext"],na.rm=T)*100,uiw=apply(Turn_prop[,,i,"Local ext"],2,sd,na.rm=T)*100, lwd=2, add=T, pch=pchV[i], col=ColV[i])
} 
legend("topright",c("Plants","Herbivores","Carnivores"),lwd=3,lty=1,col=ColV[4:6],pch=pchV[4:6], bty='n', cex=1.2)
mtext("f",3,line=0.25,cex=1.5,adj=-.21)

dev.off()
