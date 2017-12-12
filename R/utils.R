getSBMparams <- function(g, groups){
    Q <- length(unique(groups))
    l.mat <- matrix(0,Q,Q)
    pi.mat <- matrix(0,Q,Q)
    alpha.vec <- 1:Q
    adj.mat <- get.adjacency(g)
    rownames(adj.mat) <- 1:Q
    for(q in 1:Q){
      ind.q <- which(rownames(adj.mat)==q)
      for(l in 1:Q){
        ind.l <- which(rownames(adj.mat)==l)
        l.mat[q,l] <- sum(adj.mat[ind.q,ind.l])
        pi.mat[q,l] <- sum(adj.mat[ind.q,ind.l])/(length(ind.q)*length(ind.l))
      }
    }
    colnames(l.mat) <- 1:Q
    rownames(l.mat) <- 1:Q
    colnames(pi.mat) <- 1:Q
    rownames(pi.mat) <- 1:Q
    return(list(alpha=alpha.vec,
                l=l.mat,
                pi=pi.mat,
                C=sum(adj.mat)/(length(V(g))*length(V(g)))))
}

metaweb.params <- function(g.list, groups, prior.metaweb = FALSE){
    ids.list <- lapply(g.list,get.ids)   # a list of the ids (= SBM classes), within lists for each subcommunity
    metaweb.ids <- unique(unlist(ids.list))    # a list of the SBM classes in the metaweb

    Q <- length(unique(groups))
    P.mat <- matrix(0, ncol = length(g.list), nrow = Q)
    rownames(P.mat) <- 1:Q
    if(names(g.list)){
        colnames(P.mat) <- names(g.list)
    } else{
         colnames(P.mat) <- 1:length(g.list)
    }
    
    L.array <- array(0, dim=c(Q,Q,length(g.list)))
    dimnames(L.array)[[1]] <- 1:Q
    dimnames(L.array)[[2]] <- 1:Q
    dimnames(L.array)[[3]] <- colnames(P.mat)
    
    SBMparams=lapply(g.list,getSBMparams)
    
    alpha_list <- lapply(SBMparams, function(p) p$alpha)
    L_list <- lapply(SBMparams, function(p) p$l)
    Pi_list <- lapply(SBMparams, function(p) p$pi)
    
    if (!prior.metaweb){
      Pi.array.NA <- array(NA, dim = c(Q,Q,length(g.list)))
      dimnames(Pi.array.NA)[[1]] <- 1:Q
      dimnames(Pi.array.NA)[[2]] <- 1:Q 
      dimnames(Pi.array.NA)[[3]] <- colnames(P.mat)
      
      for(i in 1:length(g.list)){
        P.mat[names(alpha_list[[i]]),i] <- alpha_list[[i]]
        L.array[rownames(L_list[[i]]),colnames(L_list[[i]]),i] <- L_list[[i]]
        Pi.array.NA[rownames(Pi_list[[i]]), colnames(Pi_list[[i]]), i] <- Pi_list[[i]]
      }
      return(list(P.mat=P.mat, L.array=L.array, Pi.array=Pi.array.NA))
    }
    else {
      Pi.metaweb <- getSBMparams(g.metaweb)$pi
      Pi.array.metaweb <- array(rep(Pi.metaweb, length(g.list)), dim = c(Q,Q,length(g.list)))
      dimnames(Pi.array.metaweb)[[1]] <- metaweb.ids
      dimnames(Pi.array.metaweb)[[2]] <- metaweb.ids 
      dimnames(Pi.array.metaweb)[[3]] <- colnames(P.mat)
      
      for(i in 1:length(g.list)){
        P.mat[names(alpha_list[[i]]),i] <- alpha_list[[i]]
        L.array[rownames(L_list[[i]]),colnames(L_list[[i]]),i] <- L_list[[i]]
        Pi.array.metaweb[rownames(Pi_list[[i]]), colnames(Pi_list[[i]]), i] <- Pi_list[[i]]
      }
      return(list(P.mat=P.mat,L.array=L.array, Pi.array = Pi.array.metaweb))
      
    }  C EST QUOI CES APPELS A names / rownames / etc...???
  }
}

