# function which takes all model matrices and returns logM, logitPi and theta evaluated according to:
# Mean of the negative binomial: log(mu) = X.mu %*% a.mu + U %*% V + offset.mu
# Dispersion of the negative binomial: log(theta) = X.theta %*% a.theta + offset.theta
# Probability of the zero component: logit(P(Y=0)) = X.pi %*% alpha.pi + U %*% W + offset.pi

sim.make.matrices <- function( X.mu = NULL, a.mu=NULL , U= NULL , V = NULL , offset.mu = NULL , 
                               X.pi = NULL , alpha.pi = NULL , W = NULL , offset.pi = NULL , 
                               X.theta = NULL , a.theta = NULL, offset.theta = NULL , verbose = TRUE ) {
    
    logM <- NULL
    logitPi <- NULL
    logtheta <- NULL
    
    if ( !is.null(X.mu) & !is.null(a.mu) ){
        if ( is.null(logM) ){
            logM <- X.mu %*% a.mu
        } else {
            logM <- logM + X.mu %*% a.mu
        }
    }
    
    if ( !is.null(U) & !is.null(V) ){
        if ( is.null(logM) ){
            logM <- U %*% V
        } else {
            logM <- logM + U %*% V
        }
    }    
    
    if ( !is.null(offset.mu) ){
        if ( is.null(logM) ){
            logM <- offset.mu
        } else {
            logM <- logM + offset.mu
        }
    }

    if ( !is.null(X.pi) & !is.null(alpha.pi) ){
        if ( is.null(logitPi) ){
            logitPi <- X.pi %*% alpha.pi
        } else {
            logitPi <- logitPi + X.pi %*% alpha.pi
        }
    }
    
    if ( !is.null(U) & !is.null(W) ){
        if ( is.null(logitPi) ){
            logitPi <- U %*% W
        } else {
            logitPi <- logitPi + U %*% W
        }
    }
    
    if ( !is.null(offset.pi) ){
        if ( is.null(logitPi) ){
            logitPi <- offset.pi
        } else {
            logitPi <- logitPi + offset.pi
        }
    }
    
    if ( !is.null(X.theta) & !is.null(a.theta) ){
        if ( is.null(logtheta) ){
            logtheta <- X.theta %*% a.theta
        } else {
            logtheta <- logtheta + X.theta %*% a.theta
        }
    }
    
    if ( !is.null(offset.theta) ){
        if ( is.null(logtheta) ){
            logtheta <- offset.theta
        } else {
            logtheta <- logtheta + offset.theta
        }
    }
    
    return( list( logM = logM , logitPi = logitPi, logtheta = logtheta ) )
    
}

#' Simulate negative binomial distribution according to the following model:
#' Mean of the negative binomial: log(mu) = X.mu %*% a.mu + U %*% V + offset.mu
#' Dispersion of the negative binomial: log(theta) = X.theta %*% a.theta + offset.theta
#' Probability of the zero component: logit(P(Y=0)) = X.pi %*% alpha.pi + U %*% W + offset.pi
#' 
#' @param X.mu matrix of the model (see above, default=NULL)
#' @param X.pi matrix of the model (see above, default=NULL)
#' @param X.theta matrix of the model (see above, default=NULL)
#' @param U matrix of the model (see above, default=NULL)
#' @param V matrix of the model (see above, default=NULL)
#' @param W matrix of the model (see above, default=NULL)
#' @param offset.mu matrix of the model (see above, default=NULL)
#' @param offset.pi matrix of the model (see above, default=NULL)
#' @param offset.theta matrix of the model (see above, default=NULL)
#' @param a.mu matrix of coefficients (see above, default=NULL)
#' @param a.theta matrix of dispersion parameters (see above, default=NULL)
#' @param alpha.pi matrix of coefficients (see above, default=NULL)
#' @param sizefactors is a vector of size factors of length equal to the number of cells
#' @export

simulateNB <- function( X.mu = NULL, a.mu=NULL , U = NULL , V = NULL , offset.mu = NULL , 
                        X.pi = NULL , alpha.pi = NULL , W = NULL , offset.pi = NULL , 
                        X.theta = NULL , a.theta = NULL, offset.theta = NULL , sizefactors = NULL, 
                        verbose = NULL, no_cores = 1 ) {
        
    # evaluate matrices logM, logitPi and theta  
   
    model=sim.make.matrices ( X.mu = X.mu, a.mu = a.mu , U = U , V = V , offset.mu = offset.mu ,
                                  X.pi = X.pi , alpha.pi = alpha.pi , W = W , offset.pi = offset.pi ,
                                  X.theta = X.theta , a.theta = a.theta , offset.theta = offset.theta )
    
    inverse.logit <- binomial()$linkinv
    
    mean.nb <- as.vector( exp ( model$logM ) )    
    proba.of.zero <- as.vector ( inverse.logit ( model$logitPi ) )
    theta <- as.vector ( exp ( model$logtheta ) )
    ncells <- nrow(model$logM)
    nelements <- length(mean.nb)
    
    # simulate negative binomial with the mean matrix given by exp ( model$logM ) and dispersion
    # parameters given by exp ( model$logtheta )
    
    data.nb <- matrix(unlist( parallel::mclapply(seq(nelements), function(i) { rnbinom(1, mu = mean.nb[i] , size = theta[i]) }, mc.cores = no_cores) ), nrow = ncells )
    
    # simulate the binary dropout matrix with logit^{-1} ( logitPi ) as the probabilty of dropout
    # generation of "1" means that a dropout will be observed instead of the value
    
    data.dropout <- matrix(unlist( parallel::mclapply(seq(nelements), 
                                                      function(i) { rbinom( 1 , size =1 , prob = proba.of.zero[i] ) }, mc.cores = no_cores) ), nrow = ncells )
    
    # make final data matrix by integrating dropouts 
    
    counts <- data.nb * (1 - data.dropout)
    
    if (!is.null(sizefactors)){
        counts <- counts * matrix ( rep (sizefactors , n) , nrow = ncells ) 
    }
    
    # calculate fraction of zeros in the matrix
    
    zero.fraction <- sum ( counts == 0) / nelements    
    
    return( list ( counts = counts , data.nb = data.nb , dropouts = data.dropout , zero.fraction = zero.fraction ) )
    
}

#' Function which runs k-means on the two-dimensional projection of the data  
#' Returns all outputs of kmeansruns function from fpc 
#' @param data_kD matrix with cell's projections in k-dimensional space in lines, k is the number of columns of U
#' @param cl.range range for the possible number of clusters, default : from 2 to 10 
#' @export
kmeans_projection <- function ( data_kD , cl.range = NULL){
    
    if ( !is.null(cl.range) ){
        kmeans.result <- fpc::kmeansruns ( data = data_kD , krange = cl.range , criterion = "asw" )
    } else {
        kmeans.result <- fpc::kmeansruns ( data = data_kD , criterion = "asw" )
    }
    return ( kmeans.result )
}

#' Function which calculate rand indices comparing a reference partition to the zinb partition and  
#' to other partitions which are provided 
#' @param ref reference partition
#' @param zinb.partition partition resulting from zinb
#' @param matrix which columns are partitions resulting from other methods  
#' @export
rand.measure <- function ( ref = NULL , zinb.partition = NULL , matrix.of.partitions = NULL ) {
    
    methods.numb <- ncol ( matrix.of.partitions )
    
    zinb.rand <- mclust::adjustedRandIndex ( zinb.partition , ref )
    
    partitions.rand <- sapply(seq(1:methods.numb), function (j) { mclust::adjustedRandIndex ( matrix.of.partitions[,j] , ref ) })
    
    return ( c(zinb.rand, partitions.rand) )
}

#' Function which calculate rand indices comparing a reference partition to the zinb partition and  
#' to other partitions which are provided 
#' @param ref reference partition
#' @param zinb.partition partition resulting from zinb
#' @param matrix which columns are partitions resulting from other methods  
#' @export
silhouette.measure <- function ( ref.clust , zinb.clust , add.clust ) {
    
}

# pairwise.quality <- function(true.U,pca.result,zinb.result){
#     
#     pca.pairwisedist=dist(pca.result,method="euclidean")
#     zinb.pairwisedist=dist(zinb.result,method="euclidean")
#     true.pairwisedist=dist(true.U,method="euclidean")
#     pca.cor=cor(as.vector(pca.pairwisedist),as.vector(true.pairwisedist))
#     zinb.cor=cor(as.vector(zinb.pairwisedist),as.vector(true.pairwisedist))
#     
# }
