meta_net_plot<-function(Com,Ints,trophic,interactions=T,cut_value=0.75){
  Meta_com_init<-Com[,burnL,51:100]
  Meta_com_fin<-Com[,l+1,101:150]  
  
  if(interactions==F){Ints<-matrix(1,n,n)}
  diag(Ints)<-0
  if(trophic==T){
    colnames(Ints)<-rownames(Ints)<-c(paste("p",1:nprey),paste('h',1:npred1),paste("c",1:npred2))} else {
      colnames(Ints)<-rownames(Ints)<-paste("s",1:n)}
  
  nets_pre<-apply(Meta_com_init,2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    Int_strength_cut<-quantile(Int_strength[Int_strength>0],cut_value)#mean(Int_strength[Int_strength>0])
    Int_strength[Int_strength<Int_strength_cut]<-0
    Ints2<-1*Int_strength>0
    hold.df<-t(data.frame(Ints2[x>0,x>0]))
    net1<-graph.adjacency(hold.df)
    return(net1) 
  })
  
  
  nets_post<-apply(Meta_com_fin,2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    Int_strength_cut<-quantile(Int_strength[Int_strength>0],cut_value)#mean(Int_strength[Int_strength>0])
    Int_strength[Int_strength<Int_strength_cut]<-0
    Ints2<-1*Int_strength>0
    hold.df<-t(data.frame(Ints2[x>0,x>0]))
    net1<-graph.adjacency(hold.df)
    return(net1) 
  })

  
  mWeb1<-metaweb(nets_pre)
  mWeb2<-metaweb(nets_post)
  
  fullweb<-metaweb(list(mWeb2,mWeb1))
  vcount(fullweb)
  
  lay.mat<-data.frame(x=runif(vcount(fullweb)),y=as.numeric(factor(substring(V(fullweb)$name,1,1),levels = c("p","h","c"),ordered = T)),shape=as.numeric(factor(substring(V(fullweb)$name,1,1),levels = c("p","h","c"),ordered = T)),num=as.numeric(substring(V(fullweb)$name,3,4)), order=1:vcount(fullweb))
  lay.mat<-lay.mat%>%
    arrange(num)%>%
    group_by(y)%>%
    mutate(x=seq(0,1,length=n()))%>%
    ungroup()%>%
    mutate(y=replace(y,y==2,2.25), y=replace(y,y==1,jitter(y[y==1],amount = 0.3)))
  
  lay.mat<-lay.mat[order(lay.mat$order),]
  
  network_betaplot(mWeb2,mWeb1,layout=as.matrix(lay.mat[,1:2]),na = "forestgreen",ns = "dodgerblue3",nb = "grey",vertex.label=NA,vertex.size=8,vertex.shape="circle",edge.arrow.size=0.5,rescale=F,asp=F, ylim=c(0.5,3),xlim=c(0,1))
}

local_net_plot<-function(Com_inits,Com_final,Ints,trophic,cut_value=0.75,interactions=T){
  if(interactions==F){Ints<-matrix(1,n,n)}
  diag(Ints)<-0
  if(trophic==T){
    colnames(Ints)<-rownames(Ints)<-c(paste("p",1:nprey),paste('h',1:npred1),paste("c",1:npred2))} else {
      colnames(Ints)<-rownames(Ints)<-paste("s",1:n)}
  
  nets<-apply(cbind(Com_final,Com_inits),2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    Int_strength_cut<-quantile(Int_strength[Int_strength>0],cut_value)#mean(Int_strength[Int_strength>0])
    Int_strength[Int_strength<Int_strength_cut]<-0
    Ints2<-1*Int_strength>0
    hold.df<-t(data.frame(Ints2[x>0,x>0]))
    net1<-graph.adjacency(hold.df)
    return(net1) 
  })
  
  nets_bin<-apply(cbind(Com_final,Com_inits),2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    Int_strength_cut<-quantile(Int_strength[Int_strength>0],0.5)#mean(Int_strength[Int_strength>0])
    Int_strength[Int_strength<Int_strength_cut]<-0
    return(c(Int_strength))
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
  
  web1<-nets[[1]]
  web2<-nets[[min_dist+1]]
  
  fullweb<-metaweb(list(web1,web2))
  vcount(fullweb)
  
  if(trophic==T){
  lay.mat<-data.frame(x=runif(vcount(fullweb)),y=as.numeric(factor(substring(V(fullweb)$name,1,1),levels = c("p","h","c"),ordered = T)),shape=as.numeric(factor(substring(V(fullweb)$name,1,1),levels = c("p","h","c"),ordered = T)),num=as.numeric(substring(V(fullweb)$name,3,4)), order=1:vcount(fullweb))
  lay.mat<-lay.mat%>%
    arrange(num)%>%
    group_by(y)%>%
    mutate(x=seq(0,1,length=n()))%>%
    ungroup()%>%
    mutate(y=replace(y,y==2,2.25), y=replace(y,y==1,jitter(y[y==1],amount = 0.3)))
  
  lay.mat<-lay.mat[order(lay.mat$order),]
  par(pty='s')
  network_betaplot(web1,web2,layout=as.matrix(lay.mat[,1:2]),ns = "dodgerblue3",na = "#ff7f00",nb = "grey",vertex.label=NA,vertex.size=8,vertex.shape="circle",edge.arrow.size=0.5,rescale=F,asp=F, ylim=c(0.5,3),xlim=c(-0.1,1.1))
  } else {
    lay.mat<-as.data.frame(layout_in_circle(fullweb))
    set.seed(2)
    lay.mat<-lay.mat[sample(1:vcount(fullweb),vcount(fullweb),replace=F),]
    lay.mat$name<-V(fullweb)$name
    
    nets_namesV<-V(metaweb(list(web1,web2)))$name
    
    lay.mat_sub<-lay.mat%>%
      filter(name %in% nets_namesV)
    lay.mat_sub<-lay.mat_sub[match(nets_namesV,lay.mat_sub$name),]
    par(pty='s')
    network_betaplot(web1,web2,layout=as.matrix(lay.mat_sub[,1:2]),ns = "dodgerblue3",na = "#ff7f00",nb = "grey",vertex.label=NA,vertex.size=8,vertex.shape="circle",edge.arrow.size=0.5,rescale=F,asp=F, xlim=c(-1,1), ylim=c(-1,1))
  }
  }

nest_fun2<-function(net){
    return(NODF(import.RInSp(net))$NODF)
}

Net_ind_func<-function(Com,Ints,trophic=F,cut_value=0.75){
    trophicV<-factor(c(rep("plant",nprey),rep("herbivore",npred1),rep("predator",npred2)),levels=c("plant","herbivore","predator"),ordered = T)
  Ind_hold2<-apply(Com[,51:150,],2,function(y){
    Ind_hold<-apply(y,2,function(x){
      Int_strength<-abs(Ints*rep(x,each=n))
      Int_strength[x==0,]<-0
      Int_strength_cut<-quantile(Int_strength[Int_strength>0],cut_value)#mean(Int_strength[Int_strength>0])
      Int_strength[Int_strength<Int_strength_cut]<-0
      Ints2<-1*Int_strength>0
      hold.df<-t(data.frame(Ints2[x>0,x>0]))
      hold3<-data.frame(GenInd(hold.df))
      hold3$Nestedness<-nest_fun2(hold.df)
      if(trophic==T){hold3$Trophic_levels<-sum(tapply(x>0,trophicV,sum)>0)} else{
        hold3$Trophic_levels<-1
      }
      return(hold3)
    })
    Ind_hold2<-do.call(rbind.data.frame,Ind_hold)
    return(Ind_hold2)
  })
  

  
  Net_inds<-do.call(rbind.data.frame,Ind_hold2)
  Net_inds$time<-rep(seq(2000,7000,by=50))
  Net_inds$patch<-rep(51:150,each=101)
  return(Net_inds)
}

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

B_jack_diss<-function(pm){
  with(pm, {
    (b+c)/(a+b+c)
  })}

B_jack_diss_gains<-function(pm){
  with(pm, {
    (c)/(a+b+c)
  })}

B_jack_diss_loss<-function(pm){
  with(pm, {
    (b)/(a+b+c)
  })}

int_hist_f<-function(Com,Ints,Fpatch,trophic=F,dispersal, community){
  diag(Ints)<-0
  if(trophic==T){
    colnames(Ints)<-rownames(Ints)<-c(paste("p",1:nprey),paste('h',1:npred1),paste("c",1:npred2))} else {
      colnames(Ints)<-rownames(Ints)<-paste("s",1:n)}
  
  nets_bin<-apply(cbind(Com[,Fpatch,101],Com[,,1]),2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    Int_strength_cut<-quantile(Int_strength[Int_strength>0],0.5)#mean(Int_strength[Int_strength>0])
    Int_strength[Int_strength<Int_strength_cut]<-0
    return(c(Int_strength))
  })
  link_dis<-(as.matrix(vegdist(t(nets_bin),method = "jaccard",binary = T))[,1])[-1]
  min_dist<-which(link_dis==min(link_dis))
  
  Fints<-abs(Ints*rep(Com[,Fpatch,101],each=n))
  F_label<-1:length(Fints)
  F_label<-F_label[Fints>0]
  Fints<-Fints[Fints>0]
  
  Iints<-abs(Ints*rep(Com[,min_dist[1],1],each=n))
  I_label<-1:length(Iints)
  I_label<-I_label[Iints>0]
  Iints<-Iints[Iints>0]
  
Int.df<-data.frame(Real_ints=c(Iints,Fints),Label=c(I_label,F_label),Time=c(rep("Pre-change",length(Iints)),rep("Post-change",length((Fints)))),Dispersal = paste("Dispersal =",dispersal),Community = community)
return(Int.df)           
}
