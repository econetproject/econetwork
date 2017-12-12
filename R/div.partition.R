div.partition <- function(metaweb.array,eta=2,framework=c('RLC','Tu'),type=c('P','L','Pi')){
   if(framework=='RLC'){
    if(type=='P'){
      P.mat=metaweb.array$P.mat
      P.mat=P.mat/sum(P.mat)
      meta.P=metacommunity(P.mat)
      alphas=as.matrix(norm_sub_alpha(meta.P,qs = eta)[,c(6,7)])
      alphas.vec=as.numeric(alphas[,2])
      names(alphas.vec)=alphas[,1]
      beta=as.numeric(norm_meta_beta(meta.P,qs = eta)[7])
      gamma=as.numeric(meta_gamma(meta.P,qs=eta)[7])
      
      raw.sub.beta=raw_sub_beta(meta.P,qs=eta)
      distinc.beta=unlist(raw.sub.beta[,7])
      names(distinc.beta)=as.matrix(raw.sub.beta[6])
      
      return(list(Alphas=alphas.vec,distinctiveness.beta=distinc.beta,Beta=beta,Gamma=gamma))
      
    }
    if(type=='L'){
      n.groups=nrow(metaweb.array$P.mat)
      meta.links <- aperm(metaweb.array$L.array,c(2,1,3))  
      dim(meta.links) <- c(n.groups*n.groups,ncol(metaweb.array$P.mat)) 
      colnames(meta.links) <- colnames(metaweb.array$P.mat) 
      if(sum(rowSums(meta.links)>0)==ncol(meta.links)){
        meta.links=meta.links[-which(rowSums(meta.links)==0),]
      }
      meta.links=meta.links/(sum(meta.links))
      meta.L=metacommunity(meta.links)
      alphas=as.matrix(norm_sub_alpha(meta.L,qs = eta)[,c(6,7)])
      alphas.vec=as.numeric(alphas[,2])
      names(alphas.vec)=alphas[,1]
      
      beta=as.numeric(norm_meta_beta(meta.L,qs = eta)[7])
      gamma=as.numeric(meta_gamma(meta.L,qs=eta)[7])
      
      raw.sub.beta=raw_sub_beta(meta.L,qs=eta)
      distinc.beta=unlist(raw.sub.beta[,7])
      names(distinc.beta)=as.matrix(raw.sub.beta[6])
      
      return(list(Alphas=alphas.vec,distinctiveness.beta=distinc.beta,Beta=beta,Gamma=gamma))
    }
    if(type=='Pi'){
      n.groups=nrow(metaweb.array$P.mat)
      meta.Pi <- aperm(metaweb.array$Pi.array,c(2,1,3))  
      dim(meta.Pi) <- c(n.groups*n.groups,ncol(metaweb.array$P.mat)) 
      colnames(meta.Pi) <- colnames(metaweb.array$P.mat) 
      
      if(length(c(which(is.na(rowSums(meta.Pi))),which(rowSums(meta.Pi)==0)))>0){
        meta.Pi=meta.Pi[-c(which(is.na(rowSums(meta.Pi))),which(rowSums(meta.Pi)==0)),]
      }
      
      meta.Pi=meta.Pi/sum(meta.Pi)
      meta_Pi=metacommunity(meta.Pi)
      
      alphas=as.matrix(norm_sub_alpha(meta_Pi,qs = eta)[,c(6,7)])
      alphas.vec=as.numeric(alphas[,2])
      names(alphas.vec)=alphas[,1]
      beta=as.numeric(norm_meta_beta(meta_Pi,qs = eta)[7])
      gamma=as.numeric(meta_gamma(meta_Pi,qs=eta)[7])
      
      raw.sub.beta=raw_sub_beta(meta_Pi,qs=eta)
      distinc.beta=unlist(raw.sub.beta[,7])
      names(distinc.beta)=as.matrix(raw.sub.beta[6])
      
      return(list(Alphas=alphas.vec,distinctiveness.beta=distinc.beta,Beta=beta,Gamma=gamma))
    }
  }
  if(framework=='Tu'){
    if(type=='P'){
      abg=abgDecompQ(spxp = t(metaweb.array$P.mat),q = eta)
      return(list(Alphas=abg$Alphas,Beta=abg$Beta,Gamma=abg$Gamma))
    }
    if(type=='L'){
      n.groups=nrow(metaweb.array$P.mat)
      meta.links <- aperm(metaweb.array$L.array,c(2,1,3))  
      dim(meta.links) <- c(n.groups*n.groups,ncol(metaweb.array$P.mat)) 
      colnames(meta.links) <- colnames(metaweb.array$P.mat)
      if(sum(rowSums(meta.links)>0)==ncol(meta.links)){
        meta.links=meta.links[-which(rowSums(meta.links)==0),]
      }
      abg=abgDecompQ(spxp = t(meta.links),q = eta)
      return(list(Alphas=abg$Alphas,Beta=abg$Beta,Gamma=abg$Gamma))
      
    }
    if(type=='Pi'){
      n.groups=nrow(metaweb.array$P.mat)
      meta.Pi <- aperm(metaweb.array$Pi.array,c(2,1,3))  
      dim(meta.Pi) <- c(n.groups*n.groups,ncol(metaweb.array$P.mat)) 
      colnames(meta.Pi) <- colnames(metaweb.array$P.mat) 
      if(length(c(which(is.na(rowSums(meta.Pi))),which(rowSums(meta.Pi)==0)))>0){
        meta.Pi=meta.Pi[-c(which(is.na(rowSums(meta.Pi))),which(rowSums(meta.Pi)==0)),]
      }
      abg=abgDecompQ(spxp = t(meta.Pi),q = eta)
      return(list(Alphas=abg$Alphas,Beta=abg$Beta,Gamma=abg$Gamma))
    }
  }
}
