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

betalink_compare<-function(Com_inits,Com_final,Ints,prop_links=0.5, trophic, interactions=T,plot=T){
  if(interactions==F){Ints<-matrix(1,n,n)}
  diag(Ints)<-0
  if(trophic==T){
    colnames(Ints)<-rownames(Ints)<-c(paste("p",1:nprey),paste('h',1:npred1),paste("c",1:npred2))} else {
      colnames(Ints)<-rownames(Ints)<-paste("s",1:n)}
  nets<-apply(cbind(Com_final,Com_inits),2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
    return(graph.adjacency(t(Ints2[x>0,x>0])))
  })
  nets_bin<-apply(cbind(Com_final,Com_inits),2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
    return(c(Ints2))
  })
  link_dis<-(as.matrix(vegdist(t(nets_bin),method = "jaccard",binary = T))[,1])[-1]
  min_dist<-which(link_dis==min(link_dis))
  
  
  bl_dist<-lapply(min_dist,function(x){
    return(betalink2(nets[[x+1]],nets[[1]],bf = B_jack_diss))
  })
  bl_dist<-data.frame(matrix(unlist(bl_dist),nrow=length(min_dist),byrow=T))
  names(bl_dist)<-c("S","OS","WN","ST")
  
  if(length(min_dist)>1){
    min_dist<-min_dist[order(bl_dist$OS,decreasing = T)[1]]
    bl_dist<-bl_dist[order(bl_dist$OS,decreasing = T)[1],]
  }
  if(plot==T){network_betaplot(nets[[1]],nets[[min_dist+1]])}
  
  
  gains_loss<-rbind(betalink(nets[[1]],nets[[min_dist+1]],bf = B_jack_diss_gains),betalink(nets[[1]],nets[[min_dist+1]],bf = B_jack_diss_loss))
  gains_loss<-data.frame(matrix(unlist(gains_loss),2,byrow=F)) 
  names(gains_loss)<-c("S","OS","WN","ST")
  min_turn_all<-rbind(bl_dist,gains_loss)
  min_turn_all$patch<-min_dist
  min_turn_all$part<-c("All","Gain","Loss")
  
  return(min_turn_all)
}

betalink_min<-function(Com,Ints,prop_links=0.5, trophic, interactions=T,plot=T){
  hold<-apply(Com[,l+1,101:150],2,betalink_compare,Com_inits = Com[,burnL,25:150],Ints = Ints,prop_links = 0.5, trophic=trophic,interactions=interactions,plot=plot)
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

meta_net_turn<-function(Com,Ints,trophic,prop_links=0.5,interactions=T){
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


betalink2<-function (n1, n2, bf = B01) 
{
  v1 <- igraph::V(n1)$name
  v2 <- igraph::V(n2)$name
  vs <- v1[v1 %in% v2]
  beta_S <- bf(betapart(v1, v2))
  e1 <- plyr::aaply(igraph::get.edgelist(n1), 1, function(x) stringr::str_c(x, 
                                                                            collapse = "--", paste = "_"))
  e2 <- plyr::aaply(igraph::get.edgelist(n2), 1, function(x) stringr::str_c(x, 
                                                                            collapse = "--", paste = "_"))
  beta_WN <- bf(betapart(e1, e2))
  if (length(vs) >= 2) {
    sn1 <- igraph::induced.subgraph(n1, which(igraph::V(n1)$name %in% 
                                                vs))
    sn2 <- igraph::induced.subgraph(n2, which(igraph::V(n2)$name %in% 
                                                vs))
    se1 <- plyr::aaply(igraph::get.edgelist(sn1), 1, function(x) stringr::str_c(x, 
                                                                                collapse = "--", paste = "_"))
    se2 <- plyr::aaply(igraph::get.edgelist(sn2), 1, function(x) stringr::str_c(x, 
                                                                                collapse = "--", paste = "_"))
    beta_OS <- bf(betapart(se1, se2))
    beta_ST <- beta_WN - beta_OS
  }
  else {
    beta_OS <- NaN
    beta_ST <- NaN
  }
  return(list(S = beta_S, OS = beta_OS, WN = beta_WN, ST = beta_ST))
}