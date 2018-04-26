sbm.params <- function(g, groups=NULL){ # groups[i] for V(g)[i]
  
   if(is.null(groups)){# each node forms its own group if groups is NULL
   groups=V(g)$name
   names(groups)=V(g)$name
 }
  
  alpha.vec <- table(groups) #frequencies of the groups
    if(is.null(E(g)$weight)){#get the (weighted adjacency matrix)
    adj.mat <- get.adjacency(g)
    }  else {
    adj.mat <- get.adjacency(g,attr = "weight")
    }
  l.mat<-as.matrix(aggregate.Matrix(t(as.matrix(aggregate.Matrix(adj.mat,groups,fun='sum'))),groups,fun='sum'))[names(alpha.vec),names(alpha.vec)] #aggregate the adjacency matrix at a group level
  pi.mat<-l.mat/(alpha.vec%*%t(alpha.vec)) #link probability matrix
    
    return(list(alpha=alpha.vec,
                l=l.mat,
                pi=pi.mat,
                C=sum(adj.mat)/(length(V(g))*length(V(g)))))
}

## WARNING il faut que groups soit un tableau avec des names qui sont les sommets du metaweb 
metaweb.params <- function(g.list, groups, prior.metaweb=FALSE){ #get the L,Pi arrays and P mat for a list of graph
    groups.id <- unique(groups)
    if(length(names(g.list))){
        g.id <- names(g.list)
    } else{
        g.id <- 1:length(g.list)
    }
    Q <- length(groups.id) #number of different groups in the metaweb
    P.mat <- matrix(0, nrow=Q, ncol=length(g.list))
    rownames(P.mat) <- groups.id
    colnames(P.mat) <- g.id
    
    L.array <- array(0, dim=c(Q,Q,length(g.list))) #stacked adjacency matrix at a group level
    dimnames(L.array)[[1]] <- if(Q>1)  groups.id else list(groups.id)
    dimnames(L.array)[[2]] <- if(Q>1)  groups.id else list(groups.id)
    dimnames(L.array)[[3]] <- g.id
    
    SBMparams <- lapply(g.list, function(g){
        sbm.params(g, groups[V(g)$name])
    })
    alpha_list <- lapply(SBMparams, function(p) p$alpha)
    L_list <- lapply(SBMparams, function(p) p$l)
    Pi_list <- lapply(SBMparams, function(p) p$pi)
    
    if (!prior.metaweb){
        Pi.array.NA <- array(NA, dim = c(Q,Q,length(g.list)))
        dimnames(Pi.array.NA)[[1]] <- groups.id
        dimnames(Pi.array.NA)[[2]] <- groups.id 
        dimnames(Pi.array.NA)[[3]] <- g.id
        
        for(i in 1:length(g.list)){
            P.mat[names(alpha_list[[i]]),i] <- alpha_list[[i]]
            L.array[rownames(L_list[[i]]),colnames(L_list[[i]]),i] <- as.matrix(L_list[[i]][rownames(L_list[[i]]),colnames(L_list[[i]])])
            Pi.array.NA[rownames(Pi_list[[i]]), colnames(Pi_list[[i]]), i] <- as.matrix(Pi_list[[i]][rownames(Pi_list[[i]]),colnames(Pi_list[[i]])])
        }
        return(list(P.mat=P.mat, L.array=L.array, Pi.array=Pi.array.NA))
    }
    else {
        g.metaweb <- get.metaweb(g.list)
        Pi.metaweb <- sbm.params(g.metaweb,groups)$pi
        Pi.array.metaweb <- array(rep(Pi.metaweb, length(g.list)), dim = c(Q,Q,length(g.list)))
        dimnames(Pi.array.metaweb)[[1]] <- groups.id
        dimnames(Pi.array.metaweb)[[2]] <- groups.id 
        dimnames(Pi.array.metaweb)[[3]] <- g.id
        
        for(i in 1:length(g.list)){
            P.mat[names(alpha_list[[i]]),i] <- alpha_list[[i]]
            L.array[rownames(L_list[[i]]),colnames(L_list[[i]]),i] <- L_list[[i]]
            Pi.array.metaweb[rownames(Pi_list[[i]]), colnames(Pi_list[[i]]), i] <- Pi_list[[i]]
        }
        return(list(P.mat=P.mat, L.array=L.array, Pi.array=Pi.array.metaweb))
    }
}
