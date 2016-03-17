library(vegan)
library(plotrix)
library(RColorBrewer)

source("./Functions/predict_err.r")

dispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1)
reps<-50
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
    
    
    Envelope<-X[,burnL,51:100]
    EnvelopeI<-XI[,burnL,51:100]
    Envelope3<-X3[,burnL,51:100]
    EnvelopeM<-XM[,burnL,51:100]
    Shift<-X[,length(StressV),101:150]
    ShiftI<-XI[,length(StressV),101:150]
    Shift3<-X3[,length(StressV),101:150]
    ShiftM<-XM[,length(StressV),101:150]
    
    #Turnover####
    Turnover[r,d,"NoInt","Pre"]<-mean(colSums(Envelope>0))
    Turnover[r,d,"NoInt","Post"]<-mean(colSums(Shift>0))
    Turnover[r,d,"NoInt","Ext"]<-mean(colSums(Envelope>0 & Shift==0))
    Turnover[r,d,"NoInt","Loc_ext"]<-mean(colSums(Envelope>0 & rep(rowSums(Shift)>0,50)  & Shift==0))
    Turnover[r,d,"NoInt","Novel"]<-mean(colSums(Envelope==0 & Shift>0))
    
    Turnover[r,d,"Comp","Pre"]<-mean(colSums(EnvelopeI>0))
    Turnover[r,d,"Comp","Post"]<-mean(colSums(ShiftI>0))
    Turnover[r,d,"Comp","Ext"]<-mean(colSums(EnvelopeI>0 & ShiftI==0))
    Turnover[r,d,"Comp","Loc_ext"]<-mean(colSums(EnvelopeI>0 & rep(rowSums(ShiftI)>0,50)  & ShiftI==0))
    Turnover[r,d,"Comp","Novel"]<-mean(colSums(EnvelopeI==0 & ShiftI>0))
    
    Turnover[r,d,"Mixed","Pre"]<-mean(colSums(EnvelopeM>0))
    Turnover[r,d,"Mixed","Post"]<-mean(colSums(ShiftM>0))
    Turnover[r,d,"Mixed","Ext"]<-mean(colSums(EnvelopeM>0 & ShiftM==0))
    Turnover[r,d,"Mixed","Loc_ext"]<-mean(colSums(EnvelopeM>0 & rep(rowSums(ShiftM)>0,50)  & ShiftM==0))
    Turnover[r,d,"Mixed","Novel"]<-mean(colSums(EnvelopeM==0 & ShiftM>0))
    
    Turnover[r,d,"Plants","Pre"]<-mean(colSums(Envelope3[preyV,]>0))
    Turnover[r,d,"Plants","Post"]<-mean(colSums(Shift3[preyV,]>0))
    Turnover[r,d,"Plants","Ext"]<-mean(colSums(Envelope3[preyV,]>0 & Shift3[preyV,]==0))
    Turnover[r,d,"Plants","Loc_ext"]<-mean(colSums(Envelope3[preyV,]>0 & rep(rowSums(Shift3[preyV,])>0,50)  & Shift3[preyV,]==0))
    Turnover[r,d,"Plants","Novel"]<-mean(colSums(Envelope3[preyV,]==0 & Shift3[preyV,]>0))
    
    Turnover[r,d,"Herb","Pre"]<-mean(colSums(Envelope3[pred1V,]>0))
    Turnover[r,d,"Herb","Post"]<-mean(colSums(Shift3[pred1V,]>0))
    Turnover[r,d,"Herb","Ext"]<-mean(colSums(Envelope3[pred1V,]>0 & Shift3[pred1V,]==0))
    Turnover[r,d,"Herb","Loc_ext"]<-mean(colSums(Envelope3[pred1V,]>0 & rep(rowSums(Shift3[pred1V,])>0,50)  & Shift3[pred1V,]==0))
    Turnover[r,d,"Herb","Novel"]<-mean(colSums(Envelope3[pred1V,]==0 & Shift3[pred1V,]>0))
    
    Turnover[r,d,"Pred","Pre"]<-mean(colSums(Envelope3[pred2V,]>0))
    Turnover[r,d,"Pred","Post"]<-mean(colSums(Shift3[pred2V,]>0))
    Turnover[r,d,"Pred","Ext"]<-mean(colSums(Envelope3[pred2V,]>0 & Shift3[pred2V,]==0))
    Turnover[r,d,"Pred","Loc_ext"]<-mean(colSums(Envelope3[pred2V,]>0 & rep(rowSums(Shift3[pred2V,])>0,50)  & Shift3[pred2V,]==0))
    Turnover[r,d,"Pred","Novel"]<-mean(colSums(Envelope3[pred2V,]==0 & Shift3[pred2V,]>0))
    
    
    #Interaction Strength
    Int_Stength[,r,d,"NoInt"]<-colMeans(MeanInteract[2001:5000,,1])
    Int_Stength[,r,d,"Comp"]<-colMeans(MeanInteract[2001:5000,,2])
    Int_Stength[,r,d,"Mixed"]<-colMeans(MeanInteract[2001:5000,,3])
    Int_Stength[preyV,r,d,"Plants"]<-colMeans(MeanInteract[2001:5000,preyV,4])
    Int_Stength[pred1V,r,d,"Herb"]<-colMeans(MeanInteract[2001:5000,pred1V,4])
    Int_Stength[pred2V,r,d,"Pred"]<-colMeans(MeanInteract[2001:5000,pred2V,4])
    
    #interspecific correlations
    correlations_noint<-correlations_comp<-correlations_mixed<-c(cor(t(XM[,2000,])))
    correlations_plants<-c(cor(t(X3[preyV,2000,])))
    correlations_herb<-c(cor(t(X3[pred1V,2000,])))
    correlations_pred<-c(cor(t(X3[pred2V,2000,])))
    for(i in seq(2000,5000,length=100)){
      correlations_noint<-rbind(correlations_noint,c(cor(t(X[,i,]))))
      correlations_mixed<-rbind(correlations_mixed,c(cor(t(XM[,i,]))))
      correlations_comp<-rbind(correlations_comp,c(cor(t(XI[,i,]))))
      correlations_plants<-rbind(correlations_plants,c(cor(t(X3[preyV,i,]))))
      correlations_herb<-rbind(correlations_herb,c(cor(t(X3[pred1V,i,]))))
      correlations_pred<-rbind(correlations_pred,c(cor(t(X3[pred2V,i,]))))
    }
    correlations_noint<-correlations_noint[-1,]
    correlations_mixed<-correlations_mixed[-1,]
    correlations_comp<-correlations_comp[-1,]
    correlations_plants<-correlations_plants[-1,]
    correlations_herb<-correlations_herb[-1,]
    correlations_pred<-correlations_pred[-1,]
    
    Cor_sd[,,r,d,"NoInt"]<-matrix(apply(correlations_noint,2,sd,na.rm=T),n,n)
    diag(Cor_sd[,,r,d,"NoInt"])<-NA
    Cor_sd[,,r,d,"Mixed"]<-matrix(apply(correlations_mixed,2,sd,na.rm=T),n,n)
    diag(Cor_sd[,,r,d,"Mixed"])<-NA
    Cor_sd[,,r,d,"Comp"]<-matrix(apply(correlations_comp,2,sd,na.rm=T),n,n)
    diag(Cor_sd[,,r,d,"Comp"])<-NA
    Cor_sd[preyV,preyV,r,d,"Plants"]<-matrix(apply(correlations_plants,2,sd,na.rm=T),nprey,nprey)
    diag(Cor_sd[,,r,d,"Plants"])<-NA
    Cor_sd[pred1V,pred1V,r,d,"Herb"]<-matrix(apply(correlations_herb,2,sd,na.rm=T),npred1,npred1)
    diag(Cor_sd[,,r,d,"Herb"])<-NA
    Cor_sd[pred2V,pred2V,r,d,"Pred"]<-matrix(apply(correlations_pred,2,sd,na.rm=T),npred2,npred2)
    diag(Cor_sd[,,r,d,"Pred"])<-NA
    
    #variation in colonization rate####
    speedV<-array(NA, dim=c(80,3000,6))
    for(i in 1:3000){
      speedV[,i,1]<-apply(X[,i+2000,51:150]>1,1,which.max)
      speedV[,i,2]<-apply(XI[,i+2000,51:150]>1,1,which.max)
      speedV[,i,3]<-apply(XM[,i+2000,51:150]>1,1,which.max)
      speedV[preyV,i,4]<-apply(X3[preyV,i+2000,51:150]>1,1,which.max)
      speedV[pred1V,i,5]<-apply(X3[pred1V,i+2000,51:150]>1,1,which.max)
      speedV[pred2V,i,6]<-apply(X3[pred2V,i+2000,51:150]>1,1,which.max)
    }
    speedV[speedV==1]<-NA
    speed_mean<-matrix(NA,80,6)
    for(i in 1:80){
      speed_mean[i,1]<-mean(table(speedV[i,,1])[-c(1:2,length(table(speedV[i,,1])))])
      speed_mean[i,2]<-mean(table(speedV[i,,2])[-c(1:2,length(table(speedV[i,,2])))])
      speed_mean[i,3]<-mean(table(speedV[i,,3])[-c(1:2,length(table(speedV[i,,3])))])
      speed_mean[i,4]<-mean(table(speedV[i,,4])[-c(1:2,length(table(speedV[i,,4])))]) 
      speed_mean[i,5]<-mean(table(speedV[i,,5])[-c(1:2,length(table(speedV[i,,5])))]) 
      speed_mean[i,6]<-mean(table(speedV[i,,6])[-c(1:2,length(table(speedV[i,,6])))]) 
    }
    Shift.Speed[r,d,"NoInt"]<-sd(speed_mean[,1],na.rm=T)
    Shift.Speed[r,d,"Comp"]<-sd(speed_mean[,2],na.rm=T)
    Shift.Speed[r,d,"Mixed"]<-sd(speed_mean[,3],na.rm=T)
    Shift.Speed[r,d,"Plants"]<-sd(speed_mean[,4],na.rm=T)
    Shift.Speed[r,d,"Herb"]<-sd(speed_mean[,5],na.rm=T)
    Shift.Speed[r,d,"Pred"]<-sd(speed_mean[,6],na.rm=T)
    
    
    #Bray Curtis Similarity
    BC_Dist3<-BC_Dist3_PA<-BC_DistM<-BC_DistR_PA<-BC_Dist_PA<-BC_DistI_PA<-BC_Dist<-BC_DistI<-NA
    BC_DistM_PA<-BC_Dist_prey_PA<-BC_Dist_pred1_PA<-BC_Dist_pred2_PA<-BC_Dist_prey<-BC_Dist_pred1<-BC_Dist_pred2<-NA
    for(com in 1:50){
      BC_Dist[com]<-1-vegdist(rbind(Envelope[,com],Shift[,com]), binary=F)
      BC_DistI[com]<-1-vegdist(rbind(EnvelopeI[,com],ShiftI[,com]), binary=F)
      BC_DistM[com]<-1-vegdist(rbind(EnvelopeM[,com],ShiftM[,com]), binary=F)
      BC_Dist_prey[com]<-1-vegdist(rbind(Envelope3[preyV,][,com],Shift3[preyV,][,com]), binary=F)
      BC_Dist_pred1[com]<-1-vegdist(rbind(Envelope3[pred1V,][,com],Shift3[pred1V,][,com]), binary=F)
      BC_Dist_pred2[com]<-1-vegdist(rbind(Envelope3[pred2V,][,com],Shift3[pred2V,][,com]), binary=F)
      BC_Dist_PA[com]<-1-vegdist(rbind(Envelope[,com],Shift[,com]), binary=T)
      BC_DistI_PA[com]<-1-vegdist(rbind(EnvelopeI[,com],ShiftI[,com]), binary=T)
      BC_DistM_PA[com]<-1-vegdist(rbind(EnvelopeM[,com],ShiftM[,com]), binary=T)
      BC_Dist_prey_PA[com]<-1-vegdist(rbind(Envelope3[preyV,][,com],Shift3[preyV,][,com]), binary=T)
      BC_Dist_pred1_PA[com]<-1-vegdist(rbind(Envelope3[pred1V,][,com],Shift3[pred1V,][,com]), binary=T)
      BC_Dist_pred2_PA[com]<-1-vegdist(rbind(Envelope3[pred2V,][,com],Shift3[pred2V,][,com]), binary=T)
    }
    Bray.Curt[r,d,"AB","NoInt"]<-mean(BC_Dist)
    Bray.Curt[r,d,"AB","Comp"]<-mean(BC_DistI)
    Bray.Curt[r,d,"AB","Mixed"]<-mean(BC_DistM)
    Bray.Curt[r,d,"AB","Plants"]<-mean(BC_Dist_prey)
    Bray.Curt[r,d,"AB","Herb"]<-mean(BC_Dist_pred1)
    Bray.Curt[r,d,"AB","Pred"]<-mean(BC_Dist_pred2)
    Bray.Curt[r,d,"PA","NoInt"]<-mean(BC_Dist_PA)
    Bray.Curt[r,d,"PA","Comp"]<-mean(BC_DistI_PA)
    Bray.Curt[r,d,"PA","Mixed"]<-mean(BC_DistM_PA)
    Bray.Curt[r,d,"PA","Plants"]<-mean(BC_Dist_prey_PA)
    Bray.Curt[r,d,"PA","Herb"]<-mean(BC_Dist_pred1_PA)
    Bray.Curt[r,d,"PA","Pred"]<-mean(BC_Dist_pred2_PA)
    
    #Error rates
    Error_Rates[r,d,,"NoInt"]<-unlist(predict_err(Envelope[,],Shift[,]))
    Error_Rates[r,d,,"Comp"]<-unlist(predict_err(EnvelopeI[,],ShiftI[,]))
    Error_Rates[r,d,,"Mixed"]<-unlist(predict_err(EnvelopeM[,],ShiftM[,]))
    Error_Rates[r,d,,"Plants"]<-unlist(predict_err(Envelope3[preyV,],Shift3[preyV,]))
    Error_Rates[r,d,,"Herb"]<-unlist(predict_err(Envelope3[pred1V,],Shift3[pred1V,]))
    Error_Rates[r,d,,"Pred"]<-unlist(predict_err(Envelope3[pred2V,],Shift3[pred2V,]))
    
    #Extinction missorder
    Extinct<-colSums(apply(X,1,rowSums)!=0)
    ExtinctI<-colSums(apply(XI,1,rowSums)!=0)
    Extinct3<-colSums(apply(X3,1,rowSums)!=0)
    ExtinctM<-colSums(apply(XM,1,rowSums)!=0)
    
    Extinctprey<-Extinct3[preyV]
    Extinctpred1<-Extinct3[pred1V]
    Extinctpred2<-Extinct3[pred2V]
    
    Ext_Miss[r,d,"NoInt"]<-(sum(abs(seq(1:sum(Extinct>burnL))-order(Extinct[Extinct>burnL])))/2)/sum(Extinct>burnL)
    Ext_Miss[r,d,"Comp"]<-(sum(abs(seq(1:sum(ExtinctI>burnL))-order(ExtinctI[ExtinctI>burnL])))/2)/sum(ExtinctI>burnL)
    Ext_Miss[r,d,"Mixed"]<-(sum(abs(seq(1:sum(ExtinctM>burnL))-order(ExtinctM[ExtinctM>burnL])))/2)/sum(ExtinctM>burnL)
    Ext_Miss[r,d,"Plants"]<-mean(abs(order(T_Opt3[Extinct3>burnL])-order(Extinct3[Extinct3>burnL]))[1:sum((Extinct3[preyV]>burnL))]/2)
    Ext_Miss[r,d,"Herb"]<-mean(abs(order(T_Opt3[Extinct3>burnL])-order(Extinct3[Extinct3>burnL]))[(sum((Extinct3[preyV]>burnL))+1):sum((Extinct3[1:64]>burnL))]/2)
    Ext_Miss[r,d,"Pred"]<-mean(abs(order(T_Opt3[Extinct3>burnL])-order(Extinct3[Extinct3>burnL]))[(sum((Extinct3[1:64]>burnL))+1):sum((Extinct3>burnL))]/2)
    
    Spatial.Insurance[r,d,"L.SR","NoInt"]<-mean(specnumber(t(X[,length(StressV),101:150])))/mean(specnumber(t(X[,burnL,51:100])))
    Spatial.Insurance[r,d,"L.SR","Comp"]<-mean(specnumber(t(XI[,length(StressV),101:150])))/mean(specnumber(t(XI[,burnL,51:100])))
    Spatial.Insurance[r,d,"L.SR","Mixed"]<-mean(specnumber(t(XM[,length(StressV),101:150])))/mean(specnumber(t(XM[,burnL,51:100])))
    Spatial.Insurance[r,d,"L.SR","Plants"]<-mean(specnumber(t(X3[preyV,length(StressV),101:150])))/mean(specnumber(t(X3[preyV,burnL,51:100])))
    Spatial.Insurance[r,d,"L.SR","Herb"]<-mean(specnumber(t(X3[pred1V,length(StressV),101:150])))/mean(specnumber(t(X3[pred1V,burnL,51:100])))
    Spatial.Insurance[r,d,"L.SR","Pred"]<-mean(specnumber(t(X3[pred2V,length(StressV),101:150])))/mean(specnumber(t(X3[pred2V,burnL,51:100])))
    
    Spatial.Insurance[r,d,"R.SR","NoInt"]<-specnumber(rowSums(X[,length(StressV),101:150]))/specnumber(rowSums(X[,burnL,51:100]))
    Spatial.Insurance[r,d,"R.SR","Comp"]<-specnumber(rowSums(XI[,length(StressV),101:150]))/specnumber(rowSums(XI[,burnL,51:100]))
    Spatial.Insurance[r,d,"R.SR","Mixed"]<-specnumber(rowSums(XM[,length(StressV),101:150]))/specnumber(rowSums(XM[,burnL,51:100]))
    Spatial.Insurance[r,d,"R.SR","Plants"]<-specnumber(rowSums(X3[preyV,length(StressV),101:150]))/specnumber(rowSums(X3[preyV,burnL,51:100]))
    Spatial.Insurance[r,d,"R.SR","Herb"]<-specnumber(rowSums(X3[pred1V,length(StressV),101:150]))/specnumber(rowSums(X3[pred1V,burnL,51:100]))
    Spatial.Insurance[r,d,"R.SR","Pred"]<-specnumber(rowSums(X3[pred2V,length(StressV),101:150]))/specnumber(rowSums(X3[pred2V,burnL,51:100]))
    
    Spatial.Insurance[r,d,"B.mass","NoInt"]<-mean(colSums(X[,length(StressV),101:150]))/mean(colSums(X[,burnL,51:100]))
    Spatial.Insurance[r,d,"B.mass","Comp"]<-mean(colSums(XI[,length(StressV),101:150]))/mean(colSums(XI[,burnL,51:100]))
    Spatial.Insurance[r,d,"B.mass","Mixed"]<-mean(colSums(XM[,length(StressV),101:150]))/mean(colSums(XM[,burnL,51:100]))
    Spatial.Insurance[r,d,"B.mass","Plants"]<-mean(colSums(X3[preyV,length(StressV),101:150]))/mean(colSums(X3[preyV,burnL,51:100]))
    Spatial.Insurance[r,d,"B.mass","Herb"]<-mean(colSums(X3[pred1V,length(StressV),101:150]))/mean(colSums(X3[pred1V,burnL,51:100]))
    Spatial.Insurance[r,d,"B.mass","Pred"]<-mean(colSums(X3[pred2V,length(StressV),101:150]))/mean(colSums(X3[pred2V,burnL,51:100]))
    
    Spatial.Insurance[r,d,"L.CV","NoInt"]<-mean(apply(X[,(burnL+1):(burnL+Tmax),101:150],3,function(x){sd(x)/mean(x)}))
    Spatial.Insurance[r,d,"L.CV","Comp"]<-mean(apply(XI[,(burnL+1):(burnL+Tmax),101:150],3,function(x){sd(x)/mean(x)}))
    Spatial.Insurance[r,d,"L.CV","Mixed"]<-mean(apply(XM[,(burnL+1):(burnL+Tmax),101:150],3,function(x){sd(x)/mean(x)}))
    Spatial.Insurance[r,d,"L.CV","Plants"]<-mean(apply(X3[preyV,(burnL+1):(burnL+Tmax),101:150],3,function(x){sd(x)/mean(x)}))
    Spatial.Insurance[r,d,"L.CV","Herb"]<-mean(apply(X3[pred1V,(burnL+1):(burnL+Tmax),101:150],3,function(x){sd(x)/mean(x)}))
    Spatial.Insurance[r,d,"L.CV","Pred"]<-mean(apply(X3[pred2V,(burnL+1):(burnL+Tmax),101:150],3,function(x){sd(x)/mean(x)}))
    
    Spatial.Insurance[r,d,"R.CV","NoInt"]<-sd(apply(X[,(burnL+1):(burnL+Tmax),101:150],2,sum))/mean(apply(X[,(burnL+1):(burnL+Tmax),101:150],2,sum))
    Spatial.Insurance[r,d,"R.CV","Comp"]<-sd(apply(XI[,(burnL+1):(burnL+Tmax),101:150],2,sum))/mean(apply(XI[,(burnL+1):(burnL+Tmax),101:150],2,sum))
    Spatial.Insurance[r,d,"R.CV","Mixed"]<-sd(apply(XM[,(burnL+1):(burnL+Tmax),101:150],2,sum))/mean(apply(XM[,(burnL+1):(burnL+Tmax),101:150],2,sum)) 
    Spatial.Insurance[r,d,"R.CV","Plants"]<-sd(apply(X3[preyV,(burnL+1):(burnL+Tmax),101:150],2,sum))/mean(apply(X3[preyV,(burnL+1):(burnL+Tmax),101:150],2,sum))
    Spatial.Insurance[r,d,"R.CV","Herb"]<-sd(apply(X3[pred1V,(burnL+1):(burnL+Tmax),101:150],2,sum))/mean(apply(X3[pred1V,(burnL+1):(burnL+Tmax),101:150],2,sum))
    Spatial.Insurance[r,d,"R.CV","Pred"]<-sd(apply(X3[pred2V,(burnL+1):(burnL+Tmax),101:150],2,sum))/mean(apply(X3[pred2V,(burnL+1):(burnL+Tmax),101:150],2,sum))
  }  
};close(pb)

save(dispV,Bray.Curt,Spatial.Insurance,Shift.Speed,Error_Rates,Int_Stength,Turnover,Ext_Miss,file = "Species Interactions.RData")

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
