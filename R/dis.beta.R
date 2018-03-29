dis.beta <- function(g.list,groups,eta=1,framework=c('RLC','Tu'),type=c('P','L','Pi')){
     if(is.null(names(groups))){
        stop("groups must have names (names(groups) is NULL)")
    }
  metaweb.array <- metaweb.params(g.list,groups)
  N <- ncol( metaweb.array$P.mat)
  dis <- matrix(NA, N, N)
  
  if(framework=='RLC'){
    if(type=='P'){
      spxp=metaweb.array$P.mat
      for (i in 2:N) {
        for (j in 1:(i-1)) {
          spxp.dummy <- spxp[,c(i,j)]
          spxp.dummy= spxp.dummy/sum(spxp.dummy)
          spxp.meta=metacommunity(spxp.dummy)
          dis[i, j] <- dis[j, i] <- as.numeric(norm_meta_beta(spxp.meta,qs = eta)[7])
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
      for (i in 2:N) {
        print(i)
        for (j in 1:(i-1)) {
          spxp.dummy <- spxp[,c(i,j)]
          if(sum(rowSums(spxp.dummy)>0)<nrow(spxp.dummy)){
            spxp.dummy=spxp.dummy[-which(rowSums(spxp.dummy)==0),]
          }
          spxp.dummy= spxp.dummy/sum(spxp.dummy)
          spxp.meta=metacommunity(spxp.dummy)
          dis[i, j] <- dis[j, i] <- as.numeric(norm_meta_beta(spxp.meta,qs = eta)[7])
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
      spxp=meta.Pi
      for (i in 2:N) {
        for (j in 1:(i-1)) {
          spxp.dummy <- spxp[,c(i,j)]
          
          if(length(which(is.na(rowSums( spxp.dummy))))>0){spxp.dummy=spxp.dummy[-which(is.na(rowSums( spxp.dummy))),]}
          spxp.dummy= spxp.dummy/sum(spxp.dummy)
          spxp.meta=metacommunity(spxp.dummy)
          dis[i, j] <- dis[j, i] <- as.numeric(norm_meta_beta(spxp.meta,qs = eta)[7])
        }
      }
    }
  }
  if(framework=='Tu'){
    if(type=='P'){
      spxp=t(metaweb.array$P.mat)
      for (i in 2:N) {
        for (j in 1:(i-1)) {
          spxp.dummy <- spxp[c(i,j), ]
          res <- abgDecompQ(as.matrix(spxp.dummy), q = eta, check = FALSE)$Beta
          dis[i, j] <- dis[j, i] <- res
        }
      }
    }
    if(type=='L'){
      n.groups=nrow(metaweb.array$P.mat)
      meta.links <- aperm(metaweb.array$L.array,c(2,1,3))  
      dim(meta.links) <- c(n.groups*n.groups,ncol(metaweb.array$P.mat)) 
      colnames(meta.links) <- colnames(metaweb.array$P.mat) 
      if(sum(rowSums(meta.links)>0)==ncol(meta.links)){
        meta.links=meta.links[-which(rowSums(meta.links)==0),]
      }
      spxp=t(meta.links)
      for (i in 2:N) {
        for (j in 1:(i-1)) {
          spxp.dummy <- spxp[c(i,j), ]
          res <- abgDecompQ(as.matrix(spxp.dummy), q = eta, check = FALSE)$Beta
          dis[i, j] <- dis[j, i] <- res
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
      for (i in 2:N) {
        for (j in 1:(i-1)) {
          spxp.dummy <- spxp[c(i,j), ]
          spxp.dummy=spxp.dummy[,-which(is.na(colSums(spxp.dummy)))]
          res <- abgDecompQ(as.matrix(spxp.dummy), q = eta, check = FALSE)$Beta
          dis[i, j] <- dis[j, i] <- res
        }
      }
    }
  }
  diag(dis) <- 1
  dis <- dis - 1
  row.names(dis) <- colnames(dis) <- colnames(metaweb.array$P.mat)
  return(dis)
}
