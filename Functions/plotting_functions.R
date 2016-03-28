meta_net_plot<-function(Com,Ints,trophic,prop_links=0.5,interactions=T){
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
  
  network_betaplot(mWeb2,mWeb1,layout=as.matrix(lay.mat[,1:2]),nb = "red",vertex.label=NA,vertex.size=8,vertex.shape="circle",edge.arrow.size=0.5,rescale=F,asp=F, ylim=c(0.5,3),xlim=c(0,1))
}

local_net_plot<-function(Com_inits,Com_final,Ints,trophic,prop_links=0.5,interactions=T){
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
  
  web1<-nets[[1]]
  web2<-nets[[min_dist+1]]
  
  fullweb<-metaweb(list(web1,web2))
  vcount(fullweb)
  
  lay.mat<-data.frame(x=runif(vcount(fullweb)),y=as.numeric(factor(substring(V(fullweb)$name,1,1),levels = c("p","h","c"),ordered = T)),shape=as.numeric(factor(substring(V(fullweb)$name,1,1),levels = c("p","h","c"),ordered = T)),num=as.numeric(substring(V(fullweb)$name,3,4)), order=1:vcount(fullweb))
  lay.mat<-lay.mat%>%
    arrange(num)%>%
    group_by(y)%>%
    mutate(x=seq(0,1,length=n()))%>%
    ungroup()%>%
    mutate(y=replace(y,y==2,2.25), y=replace(y,y==1,jitter(y[y==1],amount = 0.3)))
  
  lay.mat<-lay.mat[order(lay.mat$order),]
  
  network_betaplot(web1,web2,layout=as.matrix(lay.mat[,1:2]),nb = "red",vertex.label=NA,vertex.size=8,vertex.shape="circle",edge.arrow.size=0.5,rescale=F,asp=F, ylim=c(0.5,3),xlim=c(0,1))
}