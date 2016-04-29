rShift_sim<-function(){
  #functions####
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
  
  nest_fun<-function(net){
    if(ecount(net)>0){
      return(NODF(import.RInSp(get.adjacency(net,sparse=F)))$NODF)
    } else {return(0)}
  }
  
  process_nets<-function(Initial,Final,Ints,trophic=F, interactions=T,cut_value=0.5){
    if(interactions==F){Ints<-matrix(1,80,80)} 
    diag(Ints)<-0
    if(trophic==T){
      colnames(Ints)<-rownames(Ints)<-c(paste("p",1:nprey),paste('h',1:npred1),paste("c",1:npred2))} else {
        colnames(Ints)<-rownames(Ints)<-paste("s",1:n)}
    
    #generate pre warming networks
    nets_pre<-apply(Initial[,51:150],2,function(x){
      Int_strength<-abs(Ints*rep(x,each=n))
      Int_strength[x==0,]<-0
      Int_strength_cut<-quantile(Int_strength[Int_strength>0],cut_value)#mean(Int_strength[Int_strength>0])
      Int_strength[Int_strength<Int_strength_cut]<-0
      Ints2<-1*Int_strength>0
      hold.df<-t(data.frame(Ints2[x>0,x>0]))
      net1<-graph.adjacency(hold.df)
      return(net1) 
    })
    
    #generate post warming networks
    nets_post<-apply(Final[,51:150],2,function(x){
      Int_strength<-abs(Ints*rep(x,each=n))
      Int_strength[x==0,]<-0
      Int_strength_cut<-quantile(Int_strength[Int_strength>0],cut_value)#mean(Int_strength[Int_strength>0])
      Int_strength[Int_strength<Int_strength_cut]<-0
      Ints2<-1*Int_strength>0
      hold.df<-t(data.frame(Ints2[x>0,x>0]))
      net1<-graph.adjacency(hold.df)
      return(net1) 
    })
    
    #regional networks####
    regWeb_pre<-metaweb(nets_pre[51:100])
    regWeb_post<-metaweb(nets_post[51:100])
    
    #betalink
    betaLink_R<-rbind(data.frame(betalink2(regWeb_pre,regWeb_post,bf = B_jack_diss)),
                      data.frame(betalink2(regWeb_pre,regWeb_post,bf = B_jack_diss_gains)),
                      data.frame(betalink2(regWeb_pre,regWeb_post,bf = B_jack_diss_loss)))
    betaLink_R$Part<-c("All","Gain","Loss")
    
    #net indicies
    netInd_R<-data.frame(GenInd2(get.adjacency(regWeb_post,sparse = F)))/
      data.frame(GenInd2(get.adjacency(regWeb_pre,sparse = F)))
    
    netInd_R$Nestedness<-nest_fun(regWeb_post)/nest_fun(regWeb_pre)
    
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
    
    netInd_cR$Nestedness<-nest_fun(rcWeb_post)/nest_fun(rcWeb_pre)

    
    netInd_cR$Trophic_levels<-length(unique(substring(V(rcWeb_post)$name,1,1)))/length(unique(substring(V(rcWeb_pre)$name,1,1)))
    
    #Local networks
    nets_bin<-apply(cbind(Initial[,1:150],Final[,101:150]),2,function(x){
      Int_strength<-abs(Ints*rep(x,each=n))
      Int_strength[x==0,]<-0
      Int_strength_cut<-quantile(Int_strength[Int_strength>0],cut_value)#mean(Int_strength[Int_strength>0])
      Int_strength[Int_strength<Int_strength_cut]<-0
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
    
    netInd_Ln$Nestedness<-mean(unlist(sapply(51:100,function(x){nest_fun(nets_post[[x]])})/
                                            sapply(min_dist,function(x){nest_fun(nets_pre[[x]])})),na.rm=T)
    
    netInd_Ln$Trophic_levels<-mean(unlist(sapply(51:100,function(x){length(unique(substring(V(nets_post[[x]])$name,1,1)))})/
                                            sapply(min_dist,function(x){length(unique(substring(V(nets_pre[[x]])$name,1,1)))})))
    netInd_Ln<-data.frame(netInd_Ln)
    
    #local patches####
    netInd_Lp<-colMeans(unnest(data.frame(t(sapply(51:100,function(x){GenInd2(get.adjacency(nets_post[[x]],sparse = F))}))))/
                          unnest(data.frame(t(sapply(51:100,function(x){GenInd2(get.adjacency(nets_pre[[x]],sparse = F))})))),na.rm=T)
    
    
    netInd_Lp$Nestedness<-mean(unlist(sapply(51:100,function(x){nest_fun(nets_post[[x]])})/
                                        sapply(51:100,function(x){nest_fun(nets_pre[[x]])})),na.rm=T)
    
    
    netInd_Lp$Trophic_levels<-mean(unlist(sapply(51:100,function(x){length(unique(substring(V(nets_post[[x]])$name,1,1)))})/
                                            sapply(51:100,function(x){length(unique(substring(V(nets_pre[[x]])$name,1,1)))})))
    netInd_Lp<-data.frame(netInd_Lp)
    
    
    #combine results
    betaLink<-rbind(betaLink_R,betaLink_cR,betaLink_Ln)
    betaLink$Scale<-factor(rep(c("Region","Climate region","Local network"),each=3),levels=c("Region","Climate region","Local network"),ordered = T)
    betaLink$Rep<-r
    betaLink$Dispersal<-disp
    betaLink$Cut_nets<-cut_value
    
    netInd<-rbind(netInd_R,netInd_cR,netInd_Ln,netInd_Lp)
    netInd$Scale<-factor(c("Region","Climate region","Local network","Local patch"),levels=c("Region","Climate region","Local network","Local patch"),ordered = T)
    netInd$Rep<-r
    netInd$Dispersal<-disp
    netInd$Cut_nets<-cut_value
    
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
      Trophic_BL$Cut_nets<-cut_value
      
      Trophic_NI<-rbind(regNI,rcNI,LnNI)
      Trophic_NI$Rep<-r
      Trophic_NI$Dispersal<-disp
      Trophic_NI$Cut_nets<-cut_value
      
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
    speedV<-array(NA, dim=c(80,length(sampleV),4))
    if(r==1){
      for(i in 1:length(sampleV)){
        speedV[,i,1]<-apply(X_save[,51:150,i]>1,1,which.max)
        speedV[,i,2]<-apply(XI_save[,51:150,i]>1,1,which.max)
        speedV[,i,3]<-apply(XM_save[,51:150,i]>1,1,which.max)
        speedV[,i,4]<-apply(X3_save[,51:150,i]>1,1,which.max)
      }} else{
        for(i in 1:length(sampleV)){
          speedV[,i,2]<-apply(XI_save[,51:150,i]>1,1,which.max)
          speedV[,i,3]<-apply(XM_save[,51:150,i]>1,1,which.max)
          speedV[,i,4]<-apply(X3_save[,51:150,i]>1,1,which.max)
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
  
  #the simulation script####
  dispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1)
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
  trophicV<-factor(c(rep("plant",nprey),rep("herbivore",npred1),rep("predator",npred2)),levels=c("plant","herbivore","predator"),ordered = T)
  
  
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
    
    XI<-matrix(10,n,nCom)
    XM<-X3<-X<-XI
    XIt<-XMt<-X3t<-Xt<-XI
    
    sampleV<-seq(2000,5000,by=10)
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
      if(sum(l==sampleV)==1){
        samp<-which(l==sampleV)
        X_save[,,samp]<-X
        XI_save[,,samp]<-XI
        XM_save[,,samp]<-XM
        X3_save[,,samp]<-X3
      }
    }
    
    NoInt<-rbind(process_nets(Initial=Xhold,Final=X,Ints = BN,interactions = F,cut_value = 0.25),
                 process_nets(Initial=Xhold,Final=X,Ints = BN,interactions = F,cut_value = 0.5),
                 process_nets(Initial=Xhold,Final=X,Ints = BN,interactions = F,cut_value = 0.75))
    Comp<-rbind(process_nets(Initial=XIhold,Final=XI,Ints = BI,cut_value = 0.25),
                process_nets(Initial=XIhold,Final=XI,Ints = BI,cut_value = 0.5),
                process_nets(Initial=XIhold,Final=XI,Ints = BI,cut_value = 0.75))
    Mixed<-rbind(process_nets(Initial=XMhold,Final=XM,Ints = BM,cut_value = 0.25),
                 process_nets(Initial=XMhold,Final=XM,Ints = BM,cut_value = 0.5),
                 process_nets(Initial=XMhold,Final=XM,Ints = BM,cut_value = 0.75))
    FW<-rbind(process_nets(Initial=X3hold,Final=X3,Ints = B3,trophic = T,cut_value = 0.25),
              process_nets(Initial=X3hold,Final=X3,Ints = B3,trophic = T,cut_value = 0.5),
              process_nets(Initial=X3hold,Final=X3,Ints = B3,trophic = T,cut_value = 0.75))
    
    BL_temp<-rbind.data.frame(NoInt[1,1][[1]],NoInt[2,1][[1]],NoInt[3,1][[1]],
                              Comp[1,1][[1]],Comp[2,1][[1]],Comp[3,1][[1]],
                              Mixed[1,1][[1]], Mixed[2,1][[1]], Mixed[3,1][[1]],
                              FW[1,1][[1]],FW[2,1][[1]],FW[3,1][[1]])
    BL_temp$Community<-factor(rep(c("No interactions","Competitive","Mixed","Food web"), each=9*3),levels = c("No interactions","Competitive","Mixed","Food web"),ordered = T)
    NI_temp<-rbind.data.frame(NoInt[1,2][[1]],NoInt[2,2][[1]],NoInt[3,2][[1]],
                              Comp[1,2][[1]],Comp[2,2][[1]],Comp[3,2][[1]],
                              Mixed[1,2][[1]],Mixed[2,2][[1]],Mixed[3,2][[1]],
                              FW[1,2][[1]],FW[2,2][[1]],FW[3,2][[1]])
    NI_temp$Community<-factor(rep(c("No interactions","Competitive","Mixed","Food web"), each=4*3),levels = c("No interactions","Competitive","Mixed","Food web"),ordered = T)
    
    
    BL_FW_temp<-rbind.data.frame(FW[,3][[1]])
    NI_FW_temp<-rbind.data.frame(FW[,4][[1]])
    
    Shift_temp<-Shift_sd_func()
    
    
    if(d==1){
      Net_shift.df<-BL_temp
      Net_shift_troph.df<-BL_FW_temp
      Net_ind.df<-NI_temp
      Net_ind_troph.df<-NI_FW_temp
      Shift.df<-Shift_temp
    } else {
      Net_shift.df<-rbind(Net_shift.df,BL_temp)
      Net_shift_troph.df<-rbind(Net_shift_troph.df,BL_FW_temp)
      Net_ind.df<-rbind(Net_ind.df,NI_temp)
      Net_ind_troph.df<-rbind(Net_ind_troph.df,NI_FW_temp)
      Shift.df<-rbind(Shift.df,Shift_temp)
    }
  }  
  return(list(Net_shift.df,Net_shift_troph.df,Net_ind.df,Net_ind_troph.df,Shift.df))
}
