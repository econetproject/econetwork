dis.beta <- function(g.list,groups=NULL,eta=1,framework=c('RLC','Tu'),type=c('P','L','Pi')){
  
  if(sum(unlist(lapply(lapply(g.list,FUN = function(x) V(x)$name),is.null)))>0){#check if nodes have names
    stop('nodes must have names (use V(g)$name)')
  }
   if(is.null(groups)){#if groups is NULL then each node froms its own group
    groups=unique(unlist(lapply(g.list,FUN = function(g) V(g)$name)))
    names(groups)=groups
  }
  
  if(is.null(names(groups))){#check whether groups vector has nodes
    stop("groups must have names (names(groups) is NULL)")
  }
  if(prod(names(groups) %in% unique(unlist(lapply(g.list,FUN = function(g) V(g)$name))))*prod(unique(unlist(lapply(g.list,FUN = function(g) V(g)$name))) %in% names(groups))!=1){ #check if the names of groups match to the names of the metaweb
    stop("the names of groups vector do not match to the names of the metaweb")
  }

  metaweb.array <- metaweb.params(g.list,groups) #get the metaweb array
  N <- ncol( metaweb.array$P.mat)
  dis <- matrix(NA, N, N)
  
  if(framework=='RLC'){
    if(type=='P'){
      spxp=metaweb.array$P.mat 
      pb <- txtProgressBar(min = 0, max = N*(N-1)/2, style = 3) #progression bar
      comp=0
      for (i in 2:N) {
        for (j in 1:(i-1)) {
          comp=comp+1
          spxp.dummy <- spxp[,c(i,j)]
          spxp.dummy= spxp.dummy/sum(spxp.dummy)
          spxp.meta=metacommunity(spxp.dummy) #convert to class metacom for rdiversity
          dis[i, j] <- dis[j, i] <- as.numeric(norm_meta_beta(spxp.meta,qs = eta)[7])
          setTxtProgressBar(pb,comp)
        }
      }
    }
    if(type=='L'){
     n.groups=nrow(metaweb.array$P.mat)
      meta.links <- aperm(metaweb.array$L.array,c(2,1,3))  
      dim(meta.links) <- c(n.groups*n.groups,ncol(metaweb.array$P.mat)) #transform the L array in a spxp matrix
      colnames(meta.links) <- colnames(metaweb.array$P.mat) 
      if(sum(rowSums(meta.links)>0)<nrow(meta.links)){ #remove the lines that have only 0s
        meta.links=meta.links[-which(rowSums(meta.links)==0),]
      }
      spxp=meta.links
      pb <- txtProgressBar(min = 0, max = N*(N-1)/2, style = 3)#progress bar
      comp=0
      for (i in 2:N) {
        for (j in 1:(i-1)) {
          comp=comp+1
          spxp.dummy <- spxp[,c(i,j)]
          if(sum(rowSums(spxp.dummy)>0)<nrow(spxp.dummy)){#remove the lines that have only 0s (in the pairwise matrix)
            spxp.dummy=spxp.dummy[-which(rowSums(spxp.dummy)==0),]
          }
          spxp.dummy= spxp.dummy/sum(spxp.dummy)
          spxp.meta=metacommunity(spxp.dummy)
          dis[i, j] <- dis[j, i] <- as.numeric(norm_meta_beta(spxp.meta,qs = eta)[7]) #get the pairwise beta-diversity
          setTxtProgressBar(pb, comp)
        }
      }
    }
    if(type=='Pi'){
      n.groups=nrow(metaweb.array$P.mat)
      meta.Pi <- aperm(metaweb.array$Pi.array,c(2,1,3))  #transform the Pi array in a spxp matrix
      dim(meta.Pi) <- c(n.groups*n.groups,ncol(metaweb.array$P.mat)) 
      colnames(meta.Pi) <- colnames(metaweb.array$P.mat) 
      spxp=meta.Pi
      pb <- txtProgressBar(min = 0, max = N*(N-1)/2, style = 3)  #progress bar
      comp=0
      for (i in 2:N) {
        for (j in 1:(i-1)) {
          comp=comp+1
          spxp.dummy <- spxp[,c(i,j)]
          if(length(c(which(is.na(rowSums( spxp.dummy))),which(rowSums(spxp.dummy)==0)))>0){#remove the lines that have only 0s (in the pairwise matrix)
            spxp.dummy=spxp.dummy[-c(which(is.na(rowSums( spxp.dummy))),which(rowSums(spxp.dummy)==0)),]}
          spxp.dummy= spxp.dummy/sum(spxp.dummy)
          spxp.meta=metacommunity(spxp.dummy)
          dis[i, j] <- dis[j, i] <- as.numeric(norm_meta_beta(spxp.meta,qs = eta)[7])#get the pairwise beta-diversity
          setTxtProgressBar(pb, comp)
        }
      }
    }
  }
  if(framework=='Tu'){
    if(type=='P'){
      spxp=t(metaweb.array$P.mat)
      pb <- txtProgressBar(min = 0, max = N*(N-1)/2, style = 3)
      comp=0
      for (i in 2:N) {
        for (j in 1:(i-1)) {
          comp=comp+1
          spxp.dummy <- spxp[c(i,j), ]
          res <- abgDecompQ(as.matrix(spxp.dummy), q = eta, check = FALSE)$Beta
          dis[i, j] <- dis[j, i] <- res
          setTxtProgressBar(pb, comp)
        }
      }
    }
    if(type=='L'){
      n.groups=nrow(metaweb.array$P.mat)
      meta.links <- aperm(metaweb.array$L.array,c(2,1,3))  
      dim(meta.links) <- c(n.groups*n.groups,ncol(metaweb.array$P.mat)) 
      colnames(meta.links) <- colnames(metaweb.array$P.mat) 
      if(sum(rowSums(meta.links)>0)<nrow(meta.links)){
        meta.links=meta.links[-which(rowSums(meta.links)==0),]
      }
      spxp=meta.links
      pb <- txtProgressBar(min = 0, max = N*(N-1)/2, style = 3)
      comp=0
      for (i in 2:N) {
        for (j in 1:(i-1)) {
          comp=comp+1
         spxp.dummy <- spxp[,c(i,j)]
          if(sum(rowSums(spxp.dummy)>0)<nrow(spxp.dummy)){
            spxp.dummy=spxp.dummy[-which(rowSums(spxp.dummy)==0),]
          }
          spxp.dummy= spxp.dummy/sum(spxp.dummy)
          res <- abgDecompQ(as.matrix(spxp.dummy), q = eta, check = FALSE)$Beta
          dis[i, j] <- dis[j, i] <- res
          setTxtProgressBar(pb, comp)
        }
      }
    }
    if(type=='Pi'){
      n.groups=nrow(metaweb.array$P.mat)
      meta.Pi <- aperm(metaweb.array$Pi.array,c(2,1,3))  
      dim(meta.Pi) <- c(n.groups*n.groups,ncol(metaweb.array$P.mat)) 
      colnames(meta.Pi) <- colnames(metaweb.array$P.mat) 
     if(length(c(which(is.na(rowSums(meta.Pi))),which(rowSums(meta.Pi)==0)))>0){
        meta.Pi=meta.Pi[-c(which(is.na(rowSums(meta.Pi))),which(rowSums(meta.Pi)==0)),]
      }
      spxp=t(meta.Pi)
      pb <- txtProgressBar(min = 0, max = N*(N-1)/2, style = 3)
      comp=0
      for (i in 2:N) {
        for (j in 1:(i-1)) {
          comp=comp+1
        spxp.dummy <- spxp[c(i,j), ]
        if(length(c(which(is.na(rowSums( spxp.dummy))),which(rowSums(spxp.dummy)==0)))>0){
            spxp.dummy=spxp.dummy[-c(which(is.na(rowSums( spxp.dummy))),which(rowSums(spxp.dummy)==0)),]}
         spxp.dummy= spxp.dummy/sum(spxp.dummy)
         res <- abgDecompQ(as.matrix(spxp.dummy), q = eta, check = FALSE)$Beta
         dis[i, j] <- dis[j, i] <- res
         setTxtProgressBar(pb, comp)
        }
      }
    }
  }
  diag(dis) <- 1
  dis <- dis - 1
  row.names(dis) <- colnames(dis) <- colnames(metaweb.array$P.mat)
  return(dis)
}
