B_jack_diss<-function(pm){
  with(pm, {
    (b+c)/(a+b+c)
  })}

B_jack_diss_gains<-function(pm){
  with(pm, {
    (b)/(a+b+c)
  })}

B_jack_diss_loss<-function(pm){
  with(pm, {
    (c)/(a+b+c)
  })}

betalink_compare<-function(Com_inits,Com_final,Ints,prop_links=0.5, trophic, interactions=T){
  if(interactions==F){Ints<-matrix(1,n,n)}
  diag(Ints)<-0
  if(trophic==T){
  colnames(Ints)<-rownames(Ints)<-c(paste("p",1:nprey),paste('h',1:npred1),paste("c",1:npred2))} else {
    colnames(Ints)<-rownames(Ints)<-paste("s",1:n)}
  nets<-apply(cbind(Com_final,Com_inits),2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
    return(graph.adjacency(Ints2[x>0,x>0]))
  })
  bl_dist<-lapply(nets,function(x){
    return(betalink(x,nets[[1]],bf = B_jack_diss))
  })
  bl_dist<-data.frame(matrix(unlist(bl_dist),nrow=dim(Com_inits)[2]+1,byrow=T))[-1,]
  names(bl_dist)<-c("S","OS","WN","ST")
  min_dist<-which(bl_dist$WN==min(bl_dist$WN))
  if(length(min_dist)>1){
    min_dist<-min_dist[order(bl_dist[min_dist,2],decreasing = T)[1]]
  }
  network_betaplot(nets[[1]],nets[[min_dist+1]])
  min_turn<-bl_dist[min_dist,]
  
  
  gains_loss<-rbind(betalink(nets[[1]],nets[[min_dist+1]],bf = B_jack_diss_gains),betalink(nets[[1]],nets[[min_dist+1]],bf = B_jack_diss_loss))
  gains_loss<-data.frame(matrix(unlist(gains_loss),2,byrow=F)) 
  names(gains_loss)<-c("S","OS","WN","ST")
  min_turn_all<-rbind(min_turn,gains_loss)
  min_turn_all$patch<-min_dist
  min_turn_all$part<-c("All","Gain","Loss")
   
  return(min_turn_all)
}

betalink_min<-function(Com,Ints,prop_links=0.5, trophic, interactions=T){
  hold<-apply(Com[,l+1,101:150],2,betalink_compare,Com_inits = Com[,burnL,25:150],Ints = Ints,prop_links = 0.5, trophic=trophic,interactions=interactions)
  hold2<-do.call(rbind.data.frame, hold)
  names(hold2)<-c("S","OS","WN","ST","Shift","Part")
  hold2$Patch<-rep(101:150,each=3)
  hold2$Shift<-hold2$Patch-hold2$Shift-75
  means<-hold2%>%
    group_by(Part)%>%
    summarise_each(funs(mean))
  means$Patch<-NULL
  means$Scale<-"Local"
  return(means)
}

meta_net_turn<-function(Com,Ints,trophic,interactions=T){
  Meta_com_init<-Com[,burnL,51:100]
  Meta_com_fin<-Com[,l+1,101:150]  
  
  if(interactions==F){Ints<-matrix(1,n,n)}
  diag(Ints)<-0
  if(trophic==T){
    colnames(Ints)<-rownames(Ints)<-c(paste("p",1:nprey),paste('h',1:npred1),paste("c",1:npred2))} else {
      colnames(Ints)<-rownames(Ints)<-paste("s",1:n)}
  nets1<-apply(Meta_com_init,2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
    return(graph.adjacency(Ints2[x>0,x>0]))
  })
  
  nets2<-apply(Meta_com_fin,2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
    return(graph.adjacency(Ints2[x>0,x>0]))
  })
  
  mWeb1<-metaweb(nets1)
  mWeb2<-metaweb(nets2)

  network_betaplot(mWeb2,mWeb1)
  
  meta_change<-rbind(unlist(betalink(mWeb2,mWeb1,B_jack_diss)),unlist(betalink(mWeb2,mWeb1,B_jack_diss_gains)),unlist(betalink(mWeb2,mWeb1,B_jack_diss_loss)))
  meta_change<-as.data.frame(meta_change)
  meta_change$Shift<-NA   
  meta_change$Part<-c("All","Gain","Loss")
  meta_change$Scale<-"Regional"
  
  # beta1<-network_betadiversity(nets1,bf=B_jack_diss)
  # beta1_means<-colMeans(beta1[,-c(1:2)])
  # beta2<-network_betadiversity(nets2,bf=B_jack_diss)
  # beta2_means<-colMeans(beta2[,-c(1:2)])
  
  return(meta_change)}