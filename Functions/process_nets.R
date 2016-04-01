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

process_nets<-function(Com,Ints,trophic=F, interactions=T){
  if(interactions==F){Ints<-matrix(1,80,80)} 
  diag(Ints)<-0
  if(trophic==T){
    colnames(Ints)<-rownames(Ints)<-c(paste("p",1:nprey),paste('h',1:npred1),paste("c",1:npred2))} else {
      colnames(Ints)<-rownames(Ints)<-paste("s",1:n)}
  
  #generate pre warming networks
  nets_pre<-apply(Com[,burnL,51:150],2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    mean_Int_strength<-mean(Int_strength[Int_strength>0])
    Int_strength[Int_strength<mean_Int_strength]<-0
    Ints2<-1*Int_strength>0
    hold.df<-t(data.frame(Ints2[x>0,x>0]))
    net1<-graph.adjacency(hold.df)
    return(net1) 
  })
  
  #generate post warming networks
  nets_post<-apply(Com[,l+1,51:150],2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    mean_Int_strength<-mean(Int_strength[Int_strength>0])
    Int_strength[Int_strength<mean_Int_strength]<-0
    Ints2<-1*Int_strength>0
    hold.df<-t(data.frame(Ints2[x>0,x>0]))
    net1<-graph.adjacency(hold.df)
    return(net1) 
  })
  
  #regional networks####
  regWeb_pre<-metaweb(nets_pre)
  regWeb_post<-metaweb(nets_post)
  
  #betalink
  betaLink_R<-rbind(data.frame(betalink2(regWeb_pre,regWeb_post,bf = B_jack_diss)),
  data.frame(betalink2(regWeb_pre,regWeb_post,bf = B_jack_diss_gains)),
  data.frame(betalink2(regWeb_pre,regWeb_post,bf = B_jack_diss_loss)))
  betaLink_R$Part<-c("All","Gain","Loss")
  
  #net indicies
  netInd_R<-data.frame(GenInd2(get.adjacency(regWeb_post,sparse = F)))/
    data.frame(GenInd2(get.adjacency(regWeb_pre,sparse = F)))
  
  netInd_R$Trophic_levels<-length(unique(substring(V(regWeb_post)$name,1,1)))/length(unique(substring(V(regWeb_pre)$name,1,1)))
  
  #Climate region networks####
  rcWeb_pre<-metaweb(nets_pre[1:50])
  rcWeb_post<-metaweb(nets_post[51:100])
  
  #betalink
  betaLink_cR<-rbind(data.frame(betalink2(rcWeb_pre,rcWeb_post,bf = B_jack_diss)),
                     data.frame(betalink2(rcWeb_pre,rcWeb_post,bf = B_jack_diss_gains)),
                     data.frame(betalink2(rcWeb_pre,rcWeb_post,bf = B_jack_diss_loss)))
  betaLink_cR$Part<-c("All","Gain","Loss")

  #net indicies
  netInd_cR<-data.frame(GenInd2(get.adjacency(rcWeb_post,sparse = F)))/
    data.frame(GenInd2(get.adjacency(rcWeb_pre,sparse = F)))
  netInd_cR$Trophic_levels<-length(unique(substring(V(rcWeb_post)$name,1,1)))/length(unique(substring(V(rcWeb_pre)$name,1,1)))
  
  #Local networks
  nets_bin<-apply(cbind(Com[,burnL,1:150],Com[,l+1,101:150]),2,function(x){
    Int_strength<-abs(Ints*rep(x,each=n))
    Int_strength[x==0,]<-0
    mean_Int_strength<-mean(Int_strength[Int_strength>0])
    Int_strength[Int_strength<mean_Int_strength]<-0
    return(c(Int_strength))
  })
  link_dis<-(as.matrix(vegdist(t(nets_bin),method = "jaccard",binary = T))[51:150,151:200])
  min_dist<-apply(link_dis,2,function(x){min(which(x==min(x)))})
  
  #betalink
  betaLink_Ln<-rbind.data.frame(
    colMeans(data.frame(matrix(unlist(t(mapply(betalink2,nets_pre[min_dist],nets_post[51:100],MoreArgs=list(bf=B_jack_diss)))),50,4))),
    colMeans(data.frame(matrix(unlist(t(mapply(betalink2,nets_pre[min_dist],nets_post[51:100],MoreArgs=list(bf=B_jack_diss_gains)))),50,4))),
    colMeans(data.frame(matrix(unlist(t(mapply(betalink2,nets_pre[min_dist],nets_post[51:100],MoreArgs=list(bf=B_jack_diss_loss)))),50,4))))
  names(betaLink_Ln)<-c("S","OS","WN","ST")
  betaLink_Ln$Part<-c("All","Gain","Loss")
  
  #net indicies
  netInd_Ln<-colMeans(unnest(data.frame(t(sapply(51:100,function(x){GenInd2(get.adjacency(nets_post[[x]],sparse = F))}))))/
                        unnest(data.frame(t(sapply(min_dist,function(x){GenInd2(get.adjacency(nets_pre[[x]],sparse = F))})))),na.rm=T)

  netInd_Ln$Trophic_levels<-mean(unlist(sapply(51:100,function(x){length(unique(substring(V(nets_post[[x]])$name,1,1)))})/
    sapply(min_dist,function(x){length(unique(substring(V(nets_pre[[x]])$name,1,1)))})))
  netInd_Ln<-data.frame(netInd_Ln)
  
  #local patches####
  netInd_Lp<-colMeans(unnest(data.frame(t(sapply(1:100,function(x){GenInd2(get.adjacency(nets_post[[x]],sparse = F))}))))/
                        unnest(data.frame(t(sapply(1:100,function(x){GenInd2(get.adjacency(nets_pre[[x]],sparse = F))})))),na.rm=T)
  
  netInd_Lp$Trophic_levels<-mean(unlist(sapply(1:100,function(x){length(unique(substring(V(nets_post[[x]])$name,1,1)))})/
                                          sapply(1:100,function(x){length(unique(substring(V(nets_pre[[x]])$name,1,1)))})))
  netInd_Lp<-data.frame(netInd_Lp)
  

  #combine results
  betaLink<-rbind(betaLink_R,betaLink_cR,betaLink_Ln)
  betaLink$Scale<-factor(rep(c("Region","Climate region","Local network"),each=3),levels=c("Region","Climate region","Local network"),ordered = T)
  betaLink$Rep<-r
  betaLink$Dispersal<-disp
  
  netInd<-rbind(netInd_R,netInd_cR,netInd_Ln,netInd_Lp)
  netInd$Scale<-factor(c("Region","Climate region","Local network","Local patch"),levels=c("Region","Climate region","Local network","Local patch"),ordered = T)
  netInd$Rep<-r
  netInd$Dispersal<-disp
  
  #Trophic levels
  if(trophic==T){
    hold<-TN_process(pre = regWeb_pre,post = regWeb_post)
    
    regBL<-hold[[1]]
    regBL$Scale<-"Region"
    regNI<-hold[[2]]
    regNI$Scale<-"Region"
    
    hold<-TN_process(pre = rcWeb_pre,post = rcWeb_post)
    
    rcBL<-hold[[1]]
    rcBL$Scale<-"Climate region"
    rcNI<-hold[[2]]
    rcNI$Scale<-"Climate region"
    
    for(i in 1:50){
      hold<-TN_process(nets_pre[[min_dist[i]]],nets_post[[i+50]])
      hold_beta_t<-hold[[1]]
      
      hold_NI_t<-hold[[2]]
      
      if(i==1){
        hold_beta<-hold_beta_t
        hold_NI<-hold_NI_t
      } else{
        hold_beta<-rbind(hold_beta,hold_beta_t)
        hold_NI<-rbind(hold_NI,hold_NI_t)
      }
    }
    
    LnBL<-hold_beta%>%
      group_by(Trophic,Part)%>%
      summarise_each(funs(mean(.,na.rm=T)))
    LnBL$Scale<-"Local network"

    
    LnNI<-hold_NI%>%
      group_by(Trophic)%>%
      summarise_each(funs(mean(.,na.rm=T)))
    LnNI$Scale<-"Local network"
    
    Trophic_BL<-rbind(regBL,rcBL,LnBL)
    Trophic_BL$Rep<-r
    Trophic_BL$Dispersal<-disp
    Trophic_NI<-rbind(regNI,rcNI,LnNI)
    Trophic_NI$Rep<-r
    Trophic_NI$Dispersal<-disp
    
    return(list(betaLink,netInd,Trophic_BL,Trophic_NI))
    
  } else{
    return(list(betaLink,netInd))}
}

TN_process<-function(pre,post){
  #pre
  edgeTroph<-substring(ends(pre,es=1:ecount(pre)),1,1)
  
  keep<-edgeTroph[,1]=="p" & edgeTroph[,2]=="p"
  plant_pre<-subgraph.edges(pre,eids = E(pre)[keep])
  
  keep<-edgeTroph[,1]=="p" & edgeTroph[,2]=="h"
  keep<-(keep+(edgeTroph[,1]=="h" & edgeTroph[,2]=="p"))>0
  herb_pre<-subgraph.edges(pre,eids = E(pre)[keep])
  
  keep<-edgeTroph[,1]=="c" & edgeTroph[,2]=="h"
  keep<-(keep+(edgeTroph[,1]=="h" & edgeTroph[,2]=="c"))>0
  pred_pre<-subgraph.edges(pre,eids = E(pre)[keep])
  
  #post
  if(ecount(post)>0){edgeTroph<-substring(ends(post,es=1:ecount(post)),1,1)
  
  keep<-edgeTroph[,1]=="p" & edgeTroph[,2]=="p"
  plant_post<-subgraph.edges(post,eids = E(post)[keep])
  
  keep<-edgeTroph[,1]=="p" & edgeTroph[,2]=="h"
  keep<-(keep+(edgeTroph[,1]=="h" & edgeTroph[,2]=="p"))>0
  herb_post<-subgraph.edges(post,eids = E(post)[keep])
  
  keep<-edgeTroph[,1]=="c" & edgeTroph[,2]=="h"
  keep<-(keep+(edgeTroph[,1]=="h" & edgeTroph[,2]=="c"))>0
  pred_post<-subgraph.edges(post,eids = E(post)[keep])} else{
    plant_post<-post
    herb_post<-post
    pred_post<-post
  }
  
  #betalink
  beta_troph<-rbind.data.frame(betalink2(plant_pre,plant_post,bf = B_jack_diss),
                               betalink2(plant_pre,plant_post,bf = B_jack_diss_gains),
                               betalink2(plant_pre,plant_post,bf = B_jack_diss_loss),
                               betalink2(herb_pre,herb_post,bf = B_jack_diss),
                               betalink2(herb_pre,herb_post,bf = B_jack_diss_gains),
                               betalink2(herb_pre,herb_post,bf = B_jack_diss_loss),
                               betalink2(pred_pre,pred_post,bf = B_jack_diss),
                               betalink2(pred_pre,pred_post,bf = B_jack_diss_gains),
                               betalink2(pred_pre,pred_post,bf = B_jack_diss_loss))
  beta_troph$Part<-c("All","Gain","Loss")
  beta_troph$Trophic<-factor(rep(c("Plant competition","Herbivory","Predation"),each=3),levels = c("Plant competition","Herbivory","Predation"),ordered = T)
  
  netInds_trophic<-rbind(data.frame(GenInd2(get.adjacency(plant_post,sparse = F)))/data.frame(GenInd2(get.adjacency(plant_pre,sparse = F))),
                         data.frame(GenInd2(get.adjacency(herb_post,sparse = F)))/data.frame(GenInd2(get.adjacency(herb_pre,sparse = F))),
                         data.frame(GenInd2(get.adjacency(pred_post,sparse = F)))/data.frame(GenInd2(get.adjacency(pred_pre,sparse = F))))
  netInds_trophic$Trophic<-factor(c("Plant competition","Herbivory","Predation"),levels = c("Plant competition","Herbivory","Predation"),ordered = T)
  
  return(list(beta_troph,netInds_trophic))
}

GenInd2<-function (Flow = NULL, Tij = t(Flow), Import = NULL, Export = NULL, 
                   tol = 0) 
{
  if(length(Flow)==0){
    list(N = 0, T.. = 0, TST = 0, Lint = 0, 
         Ltot = 0, LD = 0, C = 0, Tijbar = 0, 
         TSTbar = 0, Cbar = 0)
  } else{
    N <- InternalNetwork(Tij, Import, Export)
    RateComp <- N$FlowToC - N$FlowFromC
    ncTij <- ncol(Tij)
    nrTij <- nrow(Tij)
    ncomp <- ncol(N$Tint)
    compNames <- rownames(N$Tint)
    intlinks <- length(which(N$Tint > tol))
    links <- length(which(Tij > tol))
    LD <- links/ncomp
    ExportSum <- sum(N$FlowTo[N$export])
    ImportSum <- sum(N$FlowFrom[N$import])
    Throughflow <- sum(N$Tint) + ImportSum - sum(RateComp[RateComp < 
                                                            0])
    Throughput <- sum(Tij)
    Avthrflow <- Throughflow/ncomp
    Connectance <- intlinks/ncomp/(ncomp - 1)
    Avlinkweight <- Throughput/links
    linkmat <- N$Tint
    linkmat[linkmat > 0] <- 1
    Cij <- matrix(nrow = ncomp, ncol = ncomp, 0)
    for (i in 1:ncomp) {
      int_i <- union(which(linkmat[i, ] > 0), which(linkmat[, 
                                                            i] > 0))
      for (j in 1:ncomp) {
        int_j <- union(which(linkmat[j, ] > 0), which(linkmat[, 
                                                              j] > 0))
        sect <- intersect(int_i, int_j)
        uni <- union(int_i, int_j)
        Cij[i, j] <- length(sect)/length(uni)
      }
    }
    Compart <- (sum(Cij) - ncomp)/ncomp/(ncomp - 1)
    list(N = ncomp, T.. = Throughput, TST = Throughflow, Lint = intlinks, 
         Ltot = links, LD = LD, C = Connectance, Tijbar = Avlinkweight, 
         TSTbar = Avthrflow, Cbar = Compart)
  }}
environment(GenInd2) <- environment(GenInd)
  

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

#variation in colonization rate####
Shift_sd_func<-function(){
  speedV<-array(NA, dim=c(80,3000,4))
  if(r==1){
    for(i in 1:3000){
      speedV[,i,1]<-apply(X[,i+2000,51:150]>1,1,which.max)
      speedV[,i,2]<-apply(XI[,i+2000,51:150]>1,1,which.max)
      speedV[,i,3]<-apply(XM[,i+2000,51:150]>1,1,which.max)
      speedV[,i,4]<-apply(X3[,i+2000,51:150]>1,1,which.max)
    }} else{
      for(i in 1:3000){
        speedV[,i,2]<-apply(XI[,i+2000,51:150]>1,1,which.max)
        speedV[,i,3]<-apply(XM[,i+2000,51:150]>1,1,which.max)
        speedV[,i,4]<-apply(X3[,i+2000,51:150]>1,1,which.max)
      }
    }
  speedV[speedV==1]<-NA
  speed_mean<-matrix(NA,80,4)
  for(i in 1:80){
    speed_mean[i,1]<-mean(table(speedV[i,,1])[-c(length(table(speedV[i,,1])))])
    speed_mean[i,2]<-mean(table(speedV[i,,2])[-c(length(table(speedV[i,,2])))])
    speed_mean[i,3]<-mean(table(speedV[i,,3])[-c(length(table(speedV[i,,3])))])
    speed_mean[i,4]<-mean(table(speedV[i,,4])[-c(length(table(speedV[i,,4])))]) 
  }
  
  Shift_sd<-data.frame(Shift_sd=c(apply(speed_mean,2,sd,na.rm=T),unlist(tapply(speed_mean[,4],trophicV,sd,na.rm=T))),
                       Interactions=c("No interactions","Competition","Mixed","Food web","Plants","Herbivores","Carnivores"),
                       Rep=r,
                       Dispersal=disp)
  return(Shift_sd)}