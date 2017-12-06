 #####function Sunyata###
#diversity index for ecological networks###
########################

#this functions depend on the following packages :
library(igraph)
library(rdiversity)

#this function takes in input a list of adjacency matrix and 
#return the adjacency matrix of the metaweb 
get.graph.list=function(adjacency.list){
  return(lapply(adjacency.list,graph_from_adjacency_matrix))
}

#this is the inverse function
get.adjacency.list=function(g.list){
  return(lapply(g.list,get.adjacency))
}


#this function build the binary metaweb from the list of adjacency matrix of local graph. It is the union of the local graphs in the sense of graph theory.
build.binary.metaweb=function(adjacency.list){
  g.list=c()
  g.metaweb= graph_from_adjacency_matrix(adjacency.list[[1]])
  for(i in 1:length(adjacency.list)){
    g_i <- graph_from_adjacency_matrix(adjacency.list[[i]]) 
    g.list <- c(g.list,list(g_i))  # iteratively adding local networks to a list
    g.metaweb <- g.metaweb %u% g_i  # iteratively add interactions to the metaweb = union of the local networks
  }
  if(!is.null(names(adjacency.list))){
  names(g.list) <- names(adjacency.list) 
  }
  adjacency.metaweb=get.adjacency(g.metaweb)
  return(adjacency.metawab=adjacency.metaweb)
}

#dummy function to get the ids vector of the nodes of a graph
get.ids=function(g){return(unique(V(g)$id))}


# This function computes the metaweb in a slighlty different meaning :
# - first, all nodes nead to have an id (which correspond to their group or 'Elton niche')
# - then, this return the community matrix (ie the frequencies of the nodes ids accross sites (matrix), 
# the frequencies of the labelled edges accross sites (array) and the probability of interactions (array). 
#Importantly, to get an estimator of the probability of interaction between two groups in a local graph, you to observe at least one node of each group
#That's why there is a missing value in the Pi matrix  if you do not observe a node of each group. Alternatively, you can put the metaweb value.

metaweb <- function(g.list, prior.metaweb = FALSE){
  
  if(sum(unlist(lapply(lapply(g.list,get.ids),is.null)))>0 | sum(unlist(lapply(lapply(g.list,get.ids),is.na)))>0){return("give ids to the nodes of each graph")}
  else{
    
    ids.list <- lapply(g.list,get.ids)   # a list of the ids (= SBM classes), within lists for each subcommunity
    metaweb.ids <- unique(unlist(ids.list))    # a list of the SBM classes in the metaweb
    P.mat <- matrix(0, ncol = length(g.list), nrow = length(metaweb.ids))
    rownames(P.mat) <- metaweb.ids
    colnames(P.mat) <- names(g.list)
    
    L.array <- array(0, dim=c(length(metaweb.ids),length(metaweb.ids),length(g.list)))
    dimnames(L.array)[[1]] <- if(length(metaweb.ids)>1)  metaweb.ids else list(metaweb.ids)
    dimnames(L.array)[[2]] <- if(length(metaweb.ids)>1)  metaweb.ids else list(metaweb.ids)
    dimnames(L.array)[[3]] <- names(g.list)
    
    SBMparms=lapply(g.list,getSBMparms)
    
    alpha_list <- lapply(SBMparms,getSBMparms_alpha)
    L_list <- lapply(SBMparms,getSBMparms_L)
    
    Pi_list <- lapply(SBMparms,getSBMparms_Pi) # list of the interaction probabilities between groups in each subcommunity
    
    if (!prior.metaweb){
      Pi.array.NA <- array(NA, dim = c(length(metaweb.ids),length(metaweb.ids),length(g.list)))
      dimnames(Pi.array.NA)[[1]] <- metaweb.ids
      dimnames(Pi.array.NA)[[2]] <- metaweb.ids 
      dimnames(Pi.array.NA)[[3]] <- names(g.list)
      
      for(i in 1:length(g.list)){
        P.mat[names(alpha_list[[i]]),i] <- alpha_list[[i]]
        L.array[rownames(L_list[[i]]),colnames(L_list[[i]]),i] <- L_list[[i]]
        Pi.array.NA[rownames(Pi_list[[i]]), colnames(Pi_list[[i]]), i] <- Pi_list[[i]]
      }
      return(list(P.mat=P.mat,L.array=L.array, Pi.array = Pi.array.NA))
    }
    else {
      Pi.metaweb <- getSBMparms_Pi(g.metaweb)
      Pi.array.metaweb <- array(rep(Pi.metaweb, length(g.list)), dim = c(length(metaweb.ids),length(metaweb.ids),length(g.list)))
      dimnames(Pi.array.metaweb)[[1]] <- metaweb.ids
      dimnames(Pi.array.metaweb)[[2]] <- metaweb.ids 
      dimnames(Pi.array.metaweb)[[3]] <- names(g.list)
      
      for(i in 1:length(g.list)){
        P.mat[names(alpha_list[[i]]),i] <- alpha_list[[i]]
        L.array[rownames(L_list[[i]]),colnames(L_list[[i]]),i] <- L_list[[i]]
        Pi.array.metaweb[rownames(Pi_list[[i]]), colnames(Pi_list[[i]]), i] <- Pi_list[[i]]
      }
      return(list(P.mat=P.mat,L.array=L.array, Pi.array = Pi.array.metaweb))
      
    }
  }
}


#good old getSBMparms function
#this function get the parameters of a SBM given the classes (id)


getSBMparms<-function(g){
  library(igraph)
  if(is.null(V(g)$id)){return("give ids to the nodes !")}
  else{
    
    ids.set=unique(V(g)$id)
    Q=length(ids.set)
    l.mat=matrix(0,Q,Q)
    pi.mat=matrix(0,Q,Q)
    
    alpha.vec=table(V(g)$id)
    adj.mat=get.adjacency(g)
    rownames(adj.mat)=V(g)$id
    
    
    for(q in 1:Q){
      ind.q=which(rownames(adj.mat)==ids.set[q])
      for(l in 1:Q){
        ind.l=which(rownames(adj.mat)==ids.set[l])
        l.mat[q,l]=sum(adj.mat[ind.q,ind.l])
        pi.mat[q,l]=sum(adj.mat[ind.q,ind.l])/(length(ind.q)*length(ind.l))
      }
    }
    
    colnames(l.mat)= ids.set
    rownames(l.mat)= ids.set
    colnames(pi.mat)= ids.set
    rownames(pi.mat)= ids.set
    
    return(list(alpha=alpha.vec,l=l.mat,pi=pi.mat,C=sum(adj.mat)/(length(V(g))*length(V(g)))))
  }
}


#dummy functions
getSBMparms_alpha=function(SBMparms){return(SBMparms$alpha)}
getSBMparms_L=function(SBMparms){return(SBMparms$l)}
getSBMparms_Pi=function(SBMparms){return(SBMparms$pi)}

##this function gives alpha, beta and gamma diversity from a metaweb.array, on P (abundances of the group), L (abundances of the links), and Pi (probability of interactions)
#Importantly, you can choose whether you use Reeve Leinster Cobbold framework or Tuomisto framework. This uses a function codes by Loïc Chalmandrier.

ABG_decomp_eta=function(metaweb.array,eta=2,framework=c('RLC','Tu'),type=c('P','L','Pi')){
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



#this function calculate a pairwise beta, understood as a dissimilarity matrix (btwn 0 and 1) between a set of graphs
Beta_dis_eta=function(metaweb.array,eta=2,framework=c('RLC','Tu'),type=c('P','L','Pi')){
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
      if(sum(rowSums(meta.links)>0)==ncol(meta.links)){
        meta.links=meta.links[-which(rowSums(meta.links)==0),]
      }
      spxp=meta.links
      for (i in 2:N) {
        for (j in 1:(i-1)) {
          spxp.dummy <- spxp[,c(i,j)]
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



######################################################################################################################################################
#                                                                                                                                                    #
#     Alpha, beta, gamma decomposition with parametrization of the dominance (q) effect                                                              #
#                                                                                                                                                    #    
#     Refs : Chalmandrier, L., Muenkemueller, T., Lavergne, S., & Thuiller, W. (2014). Ecology, 96(1), 2015, pp. 143-153                               #
#            Chao, A., Chiu, C. H., & Jost, L. (2010). Philosophical Transactions of the Royal Society B: Biological Sciences, 365(1558), 3599-3609. #
#            Leinster, T., & Cobbold, C. A. (2012). Ecology, 93(3), 477-489.                                                                         #
#                                                                                                                                                    #
######################################################################################################################################################

# Functions
## divLeinster calculates the diversity of each site of a site by species matrix according to the q parameter according to Leinster & Cobbold 2012.
## abgDecompQ performs a alpha, beta, gamma multiplicative decomposition using Leinster's diversity indices. 

## BetaDisQ calculates the pairwise beta-diversity (minus 1) between sites of a site by species matrix according to the q parameter using the afformentionned functions
##Allows a parametrization of the dominance effect

## chaoObjects is a data preparation function. It returns adequate arguments for abgDecompQ, BetaDisQ and divLeinster to perform a diversity analysis using Chao's diversity index.
### Warning: this formula only works with an ultrametric tree!

# Arguments
## spxp : sites (row) by species (cols) matrix with or without rownames and colnames.
## Z : similarity matrix used into all functions.
## check : arguments specifying if the arguments should be checked.

divLeinster <- function(spxp, Z=NULL, q=2, check = TRUE){
  #Calcul the diversity of each site of sites by species matrix. 
  #spxp columns and Z rows and columns are assumed to be in the same order.
  if (is.null(Z)) Z <- diag(ncol(spxp))
  if (check){
    if (!inherits(spxp, "matrix")) {
      stop("object \"spxp\" is not of class \"matrix\"")}  
    if (!inherits(Z, "matrix")) {
      stop("object \"Z\" is not of class \"matrix\"")}  
    if (!all(c(ncol(Z), nrow(Z)) == ncol(spxp))){
      stop("object \"Z\" and object \"spxp\" does not have matching dimensions")}
  }  
  spxp <- sweep(spxp, 1, rowSums(spxp), "/")
  Zp <- Z %*% t(spxp)
  
  if (q != 1 & q != Inf){
    mat <- t(spxp) * (Zp)^(q-1)
    mat[is.na(mat)] <- 0
    D <- colSums(mat) ^ (1/(1-q))
  }
  if (q==Inf)  {
    D <- 1/ apply(Zp, 2, max)
  }  
  if (q == 1){
    D <- apply(Zp^t(spxp), 2, function(x) 1/prod(x))
  }
  return(D)
}

abgDecompQ <- function(spxp, Z=NULL, q=2, check=TRUE) {
  #Calcul the diversity of each site of sites by species matrix. 
  #spxp columns and Z rows/cols are assumed to be in the same order.
  if (is.null(Z)) Z <- diag(ncol(spxp))
  if (check){
    if (!inherits(spxp, "matrix")) {
      stop("object \"spxp\" is not of class \"matrix\"")}  
    if (!inherits(Z, "matrix")) {
      stop("object \"Z\" is not of class \"matrix\"")}  
    if (!all(c(ncol(Z), nrow(Z)) == ncol(spxp))){
      stop("object \"Z\" and object \"spxp\" does not have matching dimensions")}
  }
  
  site.weight <- rep(1/nrow(spxp), nrow(spxp))
  spxp <- sweep(spxp, 1, rowSums(spxp), "/")
  
  gamma.ab <- colSums(sweep(spxp, 1, site.weight, "*"))
  
  Gamma <- divLeinster(t(as.matrix(gamma.ab)), Z=Z , q=q, check = FALSE)
  Alphas <- divLeinster(spxp, Z=Z , q=q, check = FALSE) 
  
  if (q != 1 & q != Inf) {
    mAlpha <- (sum(site.weight * (Alphas ^ (1 - q))))^(1 / (1 - q))
  }
  if (q==1){
    mAlpha <- exp(sum(site.weight * log(Alphas)))
  }
  if (q==Inf){
    mAlpha <- min(Alphas)
  }
  Beta <- Gamma / mAlpha
  
  names(Alphas) <- row.names(spxp)
  res <- list(Gamma=Gamma, Beta=Beta, mAlpha=mAlpha, Alphas=Alphas)
  
  return(res)
}


