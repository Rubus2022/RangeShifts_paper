predict_err<-function(predict,act){
  prediction<-predict[rowSums(predict)>0,]
  actual<-act[rowSums(predict)>0,]
  
  N<-length(actual)
  a<-sum((prediction!=0) & (actual!=0))
  b<-sum((prediction!=0) & (actual==0))
  c<-sum((prediction==0) & (actual!=0))
  d<-sum((prediction==0) & (actual==0))
  
  Prevalence<-(a+c)/N
  Ov_Diag_Pow<-(b+d)/N
  Correct_class_rate<-(a+d)/N
  Sensitivity<-a/(a+c)
  Specificity<-d/(b+d)
  F_pos<-b/(b+d)
  F_neg<-c/(a+c)
  PPP<-a/(a+b)
  NPP<-d/(c+d)
  Miss_class<-(b+c)/N
  Odds_ratio<-(a*d)/(c*b)
  Kappa<-((a+d)-(((a+c)*(a+b)+(b+d)*(c+d))/N))/(N-(((a+c)*(a+b)+(b+d)*(c+d))/N))
  TSS<-Sensitivity+Specificity-1
  return(list(Prevalence=Prevalence,Ov_Diag_Pow=Ov_Diag_Pow,Correct_class_rate=Correct_class_rate,Sensitivity=Sensitivity,Specificity=Specificity,F_pos=F_pos,F_neg=F_neg,PPP=PPP, NPP=NPP,Miss_class=Miss_class,Odds_ratio=Odds_ratio,Kappa=Kappa,TSS=TSS)) 
}