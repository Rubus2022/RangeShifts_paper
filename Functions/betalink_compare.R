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
    if(keepV>0){
      Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
      return(graph.adjacency(t(Ints2[x>0,x>0])))} else{
        hold.df<-data.frame(Ints2[x>0,x>0])
        rownames(hold.df)<-colnames(Ints)[x>0]
        return(graph.adjacency(t(hold.df)))
      }
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
  
  if(mean(link_dis)==1){min_dist<-25}
  
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
  
  gains_loss<-rbind(betalink2(nets[[1]],nets[[min_dist+1]],bf = B_jack_diss_gains),betalink2(nets[[1]],nets[[min_dist+1]],bf = B_jack_diss_loss))
  gains_loss<-data.frame(matrix(unlist(gains_loss),2,byrow=F)) 
  names(gains_loss)<-c("S","OS","WN","ST")
  min_turn_all<-rbind(bl_dist,gains_loss)
  if(bl_dist$WN==1){min_dist<-NA}
  min_turn_all$patch<-min_dist
  min_turn_all$part<-c("All","Gain","Loss")
  
  return(min_turn_all)
}

betalink_min<-function(Com,Ints,prop_links=0.5, trophic, interactions=T,plot=T){
  hold<-apply(Com[,l+1,101:150],2,betalink_compare,Com_inits = Com[,burnL,25:150],Ints = Ints,prop_links = prop_links, trophic=trophic,interactions=interactions,plot=plot)
  hold2<-do.call(rbind.data.frame, hold)
  names(hold2)<-c("S","OS","WN","ST","Shift","Part")
  hold2$Patch<-rep(101:150,each=3)
  hold2$Shift<-hold2$Patch-hold2$Shift-75
  means<-hold2%>%
    group_by(Part)%>%
    summarise_each(funs(mean(.,na.rm=T)))
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
    Int_strength[x==0,]<-0
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    if(keepV>0){
      Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
      return(graph.adjacency(t(Ints2[x>0,x>0])))} else{
        hold.df<-data.frame(Ints2[x>0,x>0])
        rownames(hold.df)<-colnames(Ints)[x>0]
        return(graph.adjacency(t(hold.df)))
      }
  })
  
  nets2<-apply(Meta_com_fin,2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    if(keepV>0){
    Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
    return(graph.adjacency(t(Ints2[x>0,x>0])))} else{
      hold.df<-data.frame(Ints2[x>0,x>0])
      rownames(hold.df)<-colnames(Ints)[x>0]
      return(graph.adjacency(t(hold.df)))
    }
  })

  mWeb1<-metaweb(nets1)
  mWeb2<-metaweb(nets2)
  
  network_betaplot(mWeb2,mWeb1)
  
  meta_change<-rbind(unlist(betalink2(mWeb2,mWeb1,B_jack_diss)),unlist(betalink2(mWeb2,mWeb1,B_jack_diss_gains)),unlist(betalink2(mWeb2,mWeb1,B_jack_diss_loss)))
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

betalink_compare_trophic<-function(Com_inits,Com_final,Ints,prop_links=0.5, trophic=T, interactions=T,plot=F){
  if(interactions==F){Ints<-matrix(1,n,n)}
  diag(Ints)<-0
  if(trophic==T){
    colnames(Ints)<-rownames(Ints)<-c(paste("p",1:nprey),paste('h',1:npred1),paste("c",1:npred2))} else {
      colnames(Ints)<-rownames(Ints)<-paste("s",1:n)}
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
  
  if(mean(link_dis)==1){min_dist<-25}
  
  nets<-apply(cbind(Com_final,Com_inits),2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    if(keepV>0){
      Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
      return(graph.adjacency(t(Ints2[x>0,x>0])))} else{
        hold.df<-data.frame(Ints2[x>0,x>0])
        rownames(hold.df)<-colnames(Ints)[x>0]
        return(graph.adjacency(t(hold.df)))
      }
  })
  
  nets_pl<-apply(cbind(Com_final,Com_inits),2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    if(keepV>0){
      Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
      Ints3<-Ints2[preyV,preyV]
      return(graph.adjacency(t(Ints3[(x[preyV]>0),(x[preyV]>0)])))} else{
        Ints3<-Ints2[preyV,preyV]
        hold.df<-data.frame(Ints3[(x[preyV]>0),(x[preyV]>0)])
        rownames(hold.df)<-colnames(Ints[,preyV])[x[preyV]>0]
        return(graph.adjacency(t(hold.df)))
      }
  })
  
  nets_h<-apply(cbind(Com_final,Com_inits),2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    if(keepV>0){
      Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
      Ints3<-Ints2[,-pred2V]
      Ints3<-Ints3[-pred2V,]
      Ints3[preyV,preyV]<-0
      return(graph.adjacency(t(Ints3[x[-pred2V]>0,x[-pred2V]>0])))} else{
        Ints3<-Ints2[,-pred2V]
        Ints3<-Ints3[-pred2V,]
        Ints3[preyV,preyV]<-0
        hold.df<-data.frame(Ints3[x[-pred2V]>0,x[-pred2V]>0])
        rownames(hold.df)<-colnames(Ints[,-pred2V])[x[-pred2V]>0]
        return(graph.adjacency(t(hold.df)))
      }
  })
  
  nets_c<-apply(cbind(Com_final,Com_inits),2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    if(keepV>0){
      Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
      Ints3<-Ints2[-preyV,]
      Ints3<-Ints3[,-preyV]
      return(graph.adjacency(t(Ints3[x[-preyV]>0,x[-preyV]>0])))} else{
        Ints3<-Ints2[-preyV,]
        Ints3<-Ints3[,-preyV]
        hold.df<-data.frame(Ints3[x[-preyV]>0,x[-preyV]>0])
        rownames(hold.df)<-colnames(Ints[,-preyV])[x[-preyV]>0]
        return(graph.adjacency(t(hold.df)))
      }
  })
  
  bl_dist<-lapply(min_dist,function(x){
    return(betalink2(nets[[x+1]],nets[[1]],bf = B_jack_diss))
  })
  bl_dist<-data.frame(matrix(unlist(bl_dist),nrow=length(min_dist),byrow=T))
  names(bl_dist)<-c("S","OS","WN","ST")
  
  if(length(min_dist)>1){
    min_dist<-min_dist[order(bl_dist$OS,decreasing = T)[1]]
    bl_dist<-bl_dist[order(bl_dist$OS,decreasing = T)[1],]
  }
  
  Food_web_dis.temp<-rbind(data.frame(betalink2(nets_pl[[min_dist+1]],nets_pl[[1]],bf=B_jack_diss)),
                           data.frame(betalink2(nets_h[[min_dist+1]],nets_h[[1]],bf=B_jack_diss)),
                           data.frame(betalink2(nets_c[[min_dist+1]],nets_c[[1]],bf=B_jack_diss)),
                           data.frame(betalink2(nets_pl[[1]],nets_pl[[min_dist+1]],bf=B_jack_diss_gains)),
                           data.frame(betalink2(nets_h[[1]],nets_h[[min_dist+1]],bf=B_jack_diss_gains)),
                           data.frame(betalink2(nets_c[[1]],nets_c[[min_dist+1]],bf=B_jack_diss_gains)),
                           data.frame(betalink2(nets_pl[[1]],nets_pl[[min_dist+1]],bf=B_jack_diss_loss)),
                           data.frame(betalink2(nets_h[[1]],nets_h[[min_dist+1]],bf=B_jack_diss_loss)),
                           data.frame(betalink2(nets_c[[1]],nets_c[[min_dist+1]],bf=B_jack_diss_loss)))
  
  Food_web_dis.temp$Level<-c("Plants","Herbivores","Predators")
  Food_web_dis.temp$Part<-rep(c("All","Gain","Loss"),each=3)
  
  return(Food_web_dis.temp)
}

betalink_min_trophic<-function(Com,Ints,prop_links=0.5,trophic=T,interactions=T,plot=F){
  hold<-apply(Com[,l+1,101:150],2,betalink_compare_trophic,Com_inits = Com[,burnL,25:150],Ints = Ints,prop_links = prop_links, trophic=T,interactions=T,plot=plot)
  hold2<-do.call(rbind.data.frame, hold)
  means<-hold2%>%
    group_by(Part,Level)%>%
    summarise_each(funs(mean(.,na.rm=T)))
  means$Scale<-"Local"
  return(means)
}

meta_net_turn_trophic<-function(Com,Ints,prop_links=0.5,trophic=T,interactions=T,plot=F){
  Meta_com_init<-Com[,burnL,51:100]
  Meta_com_fin<-Com[,l+1,101:150]  
  
  if(interactions==F){Ints<-matrix(1,n,n)}
  diag(Ints)<-0
  if(trophic==T){
    colnames(Ints)<-rownames(Ints)<-c(paste("p",1:nprey),paste('h',1:npred1),paste("c",1:npred2))} else {
      colnames(Ints)<-rownames(Ints)<-paste("s",1:n)}
  
  
  nets1_pl<-apply(Meta_com_init,2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    if(keepV>0){
      Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
      Ints3<-Ints2[preyV,preyV]
      return(graph.adjacency(t(Ints3[(x[preyV]>0),(x[preyV]>0)])))} else{
        Ints3<-Ints2[preyV,preyV]
        hold.df<-data.frame(Ints3[(x[preyV]>0),(x[preyV]>0)])
        rownames(hold.df)<-colnames(Ints[,preyV])[x[preyV]>0]
        return(graph.adjacency(t(hold.df)))
      }
  })
  
  nets1_h<-apply(Meta_com_init,2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    if(keepV>0){
      Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
      Ints3<-Ints2[,-pred2V]
      Ints3<-Ints3[-pred2V,]
      Ints3[preyV,preyV]<-0
      return(graph.adjacency(t(Ints3[x[-pred2V]>0,x[-pred2V]>0])))} else{
        Ints3<-Ints2[,-pred2V]
        Ints3<-Ints3[-pred2V,]
        Ints3[preyV,preyV]<-0
        hold.df<-data.frame(Ints3[x[-pred2V]>0,x[-pred2V]>0])
        rownames(hold.df)<-colnames(Ints[,-pred2V])[x[-pred2V]>0]
        return(graph.adjacency(t(hold.df)))
      }
  })
  
  nets1_c<-apply(Meta_com_init,2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    if(keepV>0){
      Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
      Ints3<-Ints2[-preyV,]
      Ints3<-Ints3[,-preyV]
      return(graph.adjacency(t(Ints3[x[-preyV]>0,x[-preyV]>0])))} else{
        Ints3<-Ints2[-preyV,]
        Ints3<-Ints3[,-preyV]
        hold.df<-data.frame(Ints3[x[-preyV]>0,x[-preyV]>0])
        rownames(hold.df)<-colnames(Ints[,-preyV])[x[-preyV]>0]
        return(graph.adjacency(t(hold.df)))
      }
  })
  
  nets2_pl<-apply(Meta_com_fin,2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    if(keepV>0){
      Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
      Ints3<-Ints2[preyV,preyV]
      return(graph.adjacency(t(Ints3[(x[preyV]>0),(x[preyV]>0)])))} else{
        Ints3<-Ints2[preyV,preyV]
        hold.df<-data.frame(Ints3[(x[preyV]>0),(x[preyV]>0)])
        rownames(hold.df)<-colnames(Ints[,preyV])[x[preyV]>0]
        return(graph.adjacency(t(hold.df)))
      }
  })
  
  nets2_h<-apply(Meta_com_fin,2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    if(keepV>0){
      Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
      Ints3<-Ints2[,-pred2V]
      Ints3<-Ints3[-pred2V,]
      Ints3[preyV,preyV]<-0
      return(graph.adjacency(t(Ints3[x[-pred2V]>0,x[-pred2V]>0])))} else{
        Ints3<-Ints2[,-pred2V]
        Ints3<-Ints3[-pred2V,]
        Ints3[preyV,preyV]<-0
        hold.df<-data.frame(Ints3[x[-pred2V]>0,x[-pred2V]>0])
        rownames(hold.df)<-colnames(Ints[,-pred2V])[x[-pred2V]>0]
        return(graph.adjacency(t(hold.df)))
      }
  })
  
  nets2_c<-apply(Meta_com_fin,2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    keepV<-round(sum(Int_strength>0)*prop_links)
    Ints2<-matrix(0,dim(Ints),dim(Ints))
    colnames(Ints2)<-colnames(Ints)
    rownames(Ints2)<-rownames(Ints)
    if(keepV>0){
      Ints2[order(Int_strength,decreasing = T)[1:keepV]]<-1
      Ints3<-Ints2[-preyV,]
      Ints3<-Ints3[,-preyV]
      return(graph.adjacency(t(Ints3[x[-preyV]>0,x[-preyV]>0])))} else{
        Ints3<-Ints2[-preyV,]
        Ints3<-Ints3[,-preyV]
        hold.df<-data.frame(Ints3[x[-preyV]>0,x[-preyV]>0])
        rownames(hold.df)<-colnames(Ints[,-preyV])[x[-preyV]>0]
        return(graph.adjacency(t(hold.df)))
      }
  })
  
  mWeb1_pl<-metaweb(nets1_pl)
  mWeb2_pl<-metaweb(nets2_pl)
  mWeb1_h<-metaweb(nets1_h)
  mWeb2_h<-metaweb(nets2_h)
  mWeb1_c<-metaweb(nets1_c)
  mWeb2_c<-metaweb(nets2_c)
  
  if(plot==T){
    network_betaplot(mWeb2_pl,mWeb1_pl)
    network_betaplot(mWeb2_h,mWeb1_h)
    network_betaplot(mWeb2_c,mWeb1_c)}
  
  meta_change<-rbind(data.frame(betalink(mWeb1_pl,mWeb2_pl,bf=B_jack_diss)),
                     data.frame(betalink(mWeb1_h,mWeb2_h,bf=B_jack_diss)),
                     data.frame(betalink(mWeb1_c,mWeb2_c,bf=B_jack_diss)),
                     data.frame(betalink(mWeb1_pl,mWeb2_pl,bf=B_jack_diss_gains)),
                     data.frame(betalink(mWeb1_h,mWeb2_h,bf=B_jack_diss_gains)),
                     data.frame(betalink(mWeb1_c,mWeb2_c,bf=B_jack_diss_gains)),
                     data.frame(betalink(mWeb1_pl,mWeb2_pl,bf=B_jack_diss_loss)),
                     data.frame(betalink(mWeb1_h,mWeb2_h,bf=B_jack_diss_loss)),
                     data.frame(betalink(mWeb1_c,mWeb2_c,bf=B_jack_diss_loss)))
  
  
  meta_change$Level<-c("Plants","Herbivores","Predators")
  meta_change$Part<-rep(c("All","Gain","Loss"),each=3)
  meta_change$Scale<-"Regional"
  
  return(meta_change)}
