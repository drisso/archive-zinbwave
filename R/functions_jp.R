#' Log-likelihood of the zero-inflated negative binomial model
#'
#' @param Y the counts
#' @param mu the mean parameters of the negative binomial
#' @param theta the dispersion parameters of the negative binomial
#' @param logitP0 the logit of the probability of the zero component
#' @export
#' @importFrom copula log1pexp
zinb.loglik <- function(Y, mu, theta, logitP0) {

    # log-probabilities of counts under the NB model
    logPnb <- suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))

    # contribution of zero inflation
    #lognorm <- -logitP0 - copula::log1pexp(-logitP0)
    lognorm <- - copula::log1pexp(logitP0)
    
    # log-likelihood
    sum(logPnb[Y>0]) + sum(logPnb[Y==0] + copula::log1pexp(logitP0[Y==0]-logPnb[Y==0])) + sum(lognorm)
}





zinb.make.matrices <- function( alpha , X.mu = NULL, X.pi = NULL , X.theta = NULL , Y.mu = NULL , Y.pi = NULL , offset.mu = NULL , offset.pi = NULL , offset.theta = NULL ) {

    # Parse the model (not very clean but it works)
    i <- 0
    logM <- NULL
    logitPi <- NULL
    logtheta <- NULL
    dim.alpha <- rep(0,4)
    start.alpha <- rep(NA,4)

    if (!is.null(X.mu)) {
        n <- ncol(X.mu)
        logM <- X.mu %*% alpha[(i+1):(i+n)]
        dim.alpha[1] <- n
        start.alpha[1] <- i+1
        i <- i+n
    }

    if (!is.null(Y.mu)) {
        n <- ncol(Y.mu)
        if (is.null(logM)) {
            logM <- Y.mu %*% alpha[(i+1):(i+n)]
        } else {
            logM <- logM + Y.mu %*% alpha[(i+1):(i+n)]
        }
        logitPi <- Y.pi %*% alpha[(i+1):(i+n)]
        dim.alpha[2] <- n
        start.alpha[2] <- i+1
        i <- i+n
    }

    if (!is.null(X.pi)) {
        n <- ncol(X.pi)
        if (is.null(logitPi)) {
            logitPi <- X.pi %*% alpha[(i+1):(i+n)]
        } else {
            logitPi <- logitPi + X.pi %*% alpha[(i+1):(i+n)]
        }
        dim.alpha[3] <- n
        start.alpha[3] <- i+1
        i <- i+n
    }

    if (!is.null(X.theta)) {
        n <- ncol(X.theta)
        logtheta <- X.theta %*% alpha[(i+1):(i+n)]
        dim.alpha[4] <- n
        start.alpha[4] <- i+1
    }

    if (!is.null(offset.mu)) {
        if (is.null(logM)) {
            logM <- offset.mu
        } else {
            logM <- logM + offset.mu
        }
    }

    if (!is.null(offset.pi)) {
        if (is.null(logitPi)) {
            logitPi <- offset.pi
        } else {
            logitPi <- logitPi + offset.pi
        }
    }

    if (!is.null(offset.theta)) {
        if (is.null(logtheta)) {
            logtheta <- offset.theta
        } else {
            logtheta <- logtheta + offset.theta
        }
    }

    return(list( logM=logM , logtheta=logtheta , logitPi=logitPi , dim.alpha=dim.alpha , start.alpha=start.alpha))
}

zinb.initialize <- function(datamatrix, X=NULL, V=NULL, epsilon=10){
    # TODO: add K as input parameter, and correct the function to have K latent factors
    
    n <- nrow(datamatrix)
    J <- ncol(datamatrix)
    
    if (!is.null(X)){
        X=cbind(rep(1,n),X)
    } else {
        X <- matrix(1,ncol=1,nrow=n)
    }
    
    if (!is.null(V)){
        V=cbind(rep(1,J),V)
    } else {
        V <- matrix(1,ncol=1,nrow=J)
    }
    
    # initialization for the M part
    
    beta.mu <- matrix(0,nrow=ncol(X),ncol=J)
    gamma.mu <- matrix(0,nrow=ncol(V),ncol=n)
    
    tVgamma.mu <- t(V%*%gamma.mu)
    Xbeta.mu <- X%*%beta.mu
    
    eps.beta <- epsilon/length(beta.mu)
    eps.gamma <- epsilon/length(gamma.mu)
    eps.W <- epsilon/length(W)
    eps.alpha <- epsilon/length(alpha)
    
    P <- datamatrix!=0
    Z <- 1-as.numeric(P)
    logdata <- matrix(0,dim(datamatrix))
    logdata[P] <- log(datamatrix[P])
    
    for (m in 1:nb.repeat){
        
        # TODO: be careful that glmnet standardizes the variables; to be sure we do the correct thing see for example http://stats.stackexchange.com/questions/74206/ridge-regression-results-different-in-using-lm-ridge-and-glmnet
        
        beta.mu <- matrix(unlist( parallel::mclapply(seq(J), glmnet::glmnet (X[P[,j],], logdata[P[,j],j], offset = tVgamma.mu[,j],alpha=0,lambda=eps.beta,family="gaussian"))),nrow=ncol(X))
        Xbeta.mu <- X%*%beta.mu
        
        gamma.mu <- matrix(unlist( parallel::mclapply(seq(n), glmnet::glmnet (V[P[i,],], logdata[i,P[i,]], offset = Xbeta.mu[i,],alpha=0,lambda=eps.gamma,family="gaussian"))),nrow=ncol(V))
        tVgamma.mu <- t(V%*%gamma.mu)
    }
    
    logdata.rest <- logdata
    logdata.rest[P] <- logdata[P]-Xbeta[P]-tVgamma[P]
    logdata.rest[!P] <- NA
    
    data.impute <- softImpute::softImpute(logdata.rest,lambda=sqrt(eps.W*eps.alpha))
    # TODO only keep the top K factors
    W <- (eps.alpha/eps.W)^(1/4)*data.impute$u%*%sqrt(data.impute$d)
    alpha.mu <- (eps.W/eps.alpha)^(1/4)*sqrt(data.impute$d)%*%data.impute$v
    
    # initialization of the zero inflated part
    alpha.pi <- matrix(0,nrow=ncol(W),ncol=J)
    gamma.pi <- matrix(0,nrow=ncol(V),ncol=J)
    beta.pi <- matrix(0,nrow=ncol(X),ncol=J)
    
    for (m in 1:nb.repeat){
        # TODO: in this loop, we can just merge X and W and optimize over betaalpha; we then recover beta and alpha from betaalpha after this loop (no need to retrieve alpha and beta at each iteration)
        gamma.pi <- matrix(unlist( parallel::mclapply(seq(n), glmnet::glmnet (V, Z[i,], offset = Xbeta.pi[i,]+Walpha[i,],alpha=0,lambda=eps.gamma,family="binomial"))),nrow=ncol(V))
        tVgamma.pi <- t(V%*%gamma.pi)
        
        betaalpha <- matrix(unlist( parallel::mclapply(seq(J), glmnet::glmnet (cbind(X,W), Z[,j], offset = tVgamma.pi[,j],alpha=0,lambda=eps.beta,family="binomial"))),nrow=ncol(X)+ncol(W))
        beta.pi <- betaalpha[1:ncol(X),]
        alpha.pi <- betaalpha[(ncol(X)+1):(ncol(X)+ncol(W)),]
        Xbeta.pi <- X%*%beta.pi
        Walpha.pi <- W%*%alpha.pi
        
    }
    
    # initialization of dispersion parameter
    phi <- vector(0,length=J)
}


#' Log-likelihood of the zero-inflated negative binomial model for a regression model
#'
#' This function computes the log-likelihood of a vector of counts according to a zero-inflated negative binomial model parametrized as follows:
#' Mean of the negative binomial: log(mu) = X.mu %*% a.mu + Y.mu %*%* b + offset.mu
#' Dispersion of the negative binomial: log(theta) = X.theta %*% a.theta + offset.theta
#' Probability of the zero component: logit(P(Y=0)) = X.pi %*% alpha.pi + Y.pi %*% b + offset.pi
#' Note that the b vector is shared between the mean of the negative binomial and the probability of zero.
#'
#' @param alpha the vectors of parameters c(a.mu, b, a.pi, a.theta) concatenated
#' @param Y the vector of counts
#' @param X.mu matrix of the model (see above, default=NULL)
#' @param X.pi matrix of the model (see above, default=NULL)
#' @param X.theta matrix of the model (see above, default=NULL)
#' @param Y.mu matrix of the model (see above, default=NULL)
#' @param Y.pi matrix of the model (see above, default=NULL)
#' @param offset.mu matrix of the model (see above, default=NULL)
#' @param offset.pi matrix of the model (see above, default=NULL)
#' @param offset.theta matrix of the model (see above, default=NULL)
#' @param epsilon the regularization parameter (default=0)
#' @export
zinb.loglik.regression <- function( alpha , Y, X.mu = NULL, X.pi = NULL , X.theta = NULL , Y.mu = NULL , Y.pi = NULL , offset.mu = NULL , offset.pi = NULL , offset.theta = NULL , epsilon=0) {

    # Parse the model
    r <- zinb.make.matrices(alpha=alpha, X.mu=X.mu, X.pi=X.pi, X.theta=X.theta , Y.mu = Y.mu, Y.pi=Y.pi , offset.mu=offset.mu , offset.pi=offset.pi , offset.theta=offset.theta)

    # Call the log likelihood function
    z=zinb.loglik(Y, exp(r$logM), exp(r$logtheta), r$logitPi) - epsilon*sum(alpha^2)/2
#    print(z)
    z
}


#' Gradient of the log-likelihood of the zero-inflated negative binomial model for a regression model
#'
#' This function computes the gradient of the log-likelihood of a vector of counts according to a zero-inflated negative binomial model parametrized as follows:
#' Mean of the negative binomial: log(mu) = X.mu %*% a.mu + Y.mu %*%* b + offset.mu
#' Dispersion of the negative binomial: log(theta) = X.theta %*% a.theta + offset.theta
#' Probability of the zero component: logit(P(Y=0)) = X.pi %*% alpha.pi + Y.pi %*% b + offset.pi
#' Note that the b vector is shared between the mean of the negative binomial and the probability of zero.
#'
#' @param alpha the vectors of parameters c(a.mu, b, a.pi, a.theta) concatenated
#' @param Y the vector of counts
#' @param X.mu matrix of the model (see above, default=NULL)
#' @param X.pi matrix of the model (see above, default=NULL)
#' @param X.theta matrix of the model (see above, default=NULL)
#' @param Y.mu matrix of the model (see above, default=NULL)
#' @param Y.pi matrix of the model (see above, default=NULL)
#' @param offset.mu matrix of the model (see above, default=NULL)
#' @param offset.pi matrix of the model (see above, default=NULL)
#' @param offset.theta matrix of the model (see above, default=NULL)
#' @param epsilon regularization parameter (default=0)
#' @export
gradient.zinb.loglik.regression <- function( alpha , Y, X.mu = NULL, X.pi = NULL , X.theta = NULL , Y.mu = NULL , Y.pi = NULL , offset.mu = NULL , offset.pi = NULL , offset.theta = NULL , epsilon=0) {

    # Parse the model
    r <- zinb.make.matrices(alpha=alpha, X.mu=X.mu, X.pi=X.pi, X.theta=X.theta , Y.mu = Y.mu, Y.pi=Y.pi , offset.mu=offset.mu , offset.pi=offset.pi , offset.theta=offset.theta)
    theta <- exp(r$logtheta)
    mu <- exp(r$logM)
    n <- length(Y)

    # Check zeros in the count matrix
    Y0 <- Y <= 0
    Y1 <- Y > 0
    has0 <- !is.na(match(TRUE,Y0))
    has1 <- !is.na(match(TRUE,Y1))

    # Check what we need to compute, depending on the variables over which we optimize
    need.wres.mu <- r$dim.alpha[1] >0 || r$dim.alpha[2] >0
    need.wres.pi <- r$dim.alpha[2] >0 || r$dim.alpha[3] >0
    need.wres.theta <- r$dim.alpha[4] >0

    # Compute some useful quantities
    muz <- 1/(1+exp(-r$logitPi))
    clogdens0 <- dnbinom(0, size = theta[Y0], mu = mu[Y0], log = TRUE)
    # dens0 <- muz[Y0] + exp(log(1 - muz[Y0]) + clogdens0)
    # More accurate: log(1-muz) is the following
    lognorm <- -r$logitPi - copula::log1pexp(-r$logitPi)

    dens0 <- muz[Y0] + exp(lognorm[Y0] + clogdens0)

    # Compute the partial derivatives we need
    if (need.wres.mu) {
        wres_mu <- numeric(length = n)
        if (has1) {
            wres_mu[Y1] <- Y[Y1] - mu[Y1] * (Y[Y1] + theta[Y1])/(mu[Y1] + theta[Y1])
        }
        if (has0) {
            #wres_mu[Y0] <- -exp(-log(dens0) + log(1 - muz[Y0]) + clogdens0 + r$logtheta[Y0] - log(mu[Y0] + theta[Y0]) + log(mu[Y0]))
            # more accurate:
            wres_mu[Y0] <- -exp(-log(dens0) + lognorm[Y0] + clogdens0 + r$logtheta[Y0] - log(mu[Y0] + theta[Y0]) + log(mu[Y0]))
        }
    }

    if (need.wres.pi) {
        wres_pi <- numeric(length = n)
        if (has1) {
#            wres_pi[Y1] <- -1/(1 - muz[Y1]) * linkobj$mu.eta(r$logitPi)[Y1]
            # simplification for the logit link function:
            wres_pi[Y1] <- - muz[Y1]

        }
        if (has0) {
            # wres_pi[Y0] <- (linkobj$mu.eta(r$logitPi)[Y0] - exp(clogdens0) * linkobj$mu.eta(r$logitPi)[Y0])/dens0
            # simplification for the logit link function:
            wres_pi[Y0] <- (1 - exp(clogdens0)) * muz[Y0] * (1-muz[Y0]) / dens0
        }
    }

    if (need.wres.theta) {
        wres_theta <- numeric(length = n)
        if (has1) {
            wres_theta[Y1] <- theta[Y1] * (digamma(Y[Y1] + theta[Y1]) - digamma(theta[Y1]) + r$logtheta[Y1] - log(mu[Y1] + theta[Y1]) + 1 - (Y[Y1] + theta[Y1])/(mu[Y1] + theta[Y1]) )
        }
        if (has0) {
            #wres_theta[Y0] <- theta[Y0] * ( exp(-log(dens0) + log(1 - muz[Y0]) + clogdens0) * (r$logtheta[Y0] - log(mu[Y0] + theta[Y0]) + 1 - theta[Y0]/(mu[Y0] + theta[Y0])) )
            # more accurate:
            wres_theta[Y0] <- theta[Y0] * ( exp(-log(dens0) + lognorm[Y0] + clogdens0) * (r$logtheta[Y0] - log(mu[Y0] + theta[Y0]) + 1 - theta[Y0]/(mu[Y0] + theta[Y0])) )
        }
    }

    # Make gradient
    grad <- numeric(0)
    if (r$dim.alpha[1] >0) { grad <- c(grad , colSums(wres_mu * X.mu) - epsilon*alpha[r$start.alpha[1]:(r$start.alpha[1]+r$dim.alpha[1]-1)])}
    if (r$dim.alpha[2] >0) { grad <- c(grad , colSums(wres_mu * Y.mu) + colSums(wres_pi * Y.pi) - epsilon*alpha[r$start.alpha[2]:(r$start.alpha[2]+r$dim.alpha[2]-1)] )}
    if (r$dim.alpha[3] >0) { grad <- c(grad , colSums(wres_pi * X.pi) - epsilon*alpha[r$start.alpha[3]:(r$start.alpha[3]+r$dim.alpha[3]-1)] )}
    if (r$dim.alpha[4] >0) { grad <- c(grad , colSums(wres_theta * X.theta) - epsilon*alpha[r$start.alpha[4]:(r$start.alpha[4]+r$dim.alpha[4]-1)] )}

    grad
}


#' Orthogonalize U, V, W
#'
#' Given U, V, W, find U2, V2, W2 such that UV'=U2V2' and
#' UW'=U2W2' and such that ||U2||^2+||V2||^2+||W||^2 is minimal.
#' @param U left matrix
#' @param V first right matrix
#' @param W second right matrix
#' @export
orthogonalizeJointly <- function(U, V, W) {

    # do QR of U
    U.qr <- qr (U)
    U.Q <- qr.Q (U.qr)
    U.R <- qr.R (U.qr)

    #do QR of [V;W]
    VW.qr <- qr (rbind(V,W))
    VW.Q <- qr.Q (VW.qr)
    VW.R <- qr.R (VW.qr)

    # do SVD of the U.R %*% t(V.R) matrix to have orthog %*% diag %*% orthog
    A <- svd( U.R %*% t(VW.R) )

    # orthogonalized U
    U2 <- U.Q %*% A$u %*% sqrt(diag(A$d))

    # orthogonalized V and W
    p <- nrow(V)
    VW <- VW.Q %*% A$v %*% sqrt(diag(A$d))
    V2 <- VW[1:p,]
    W2 <- VW[(p+1):(2*p),]

    list(U=U2, V=V2, W=W2)
}

#' Estimation of latent factors from a count matrix (slow version)
#' @param datamatrix the count data matrix (cells in rows, genes in columns)
#' @param k number of latent factors (default 2)
#' @param alt.number maximum number of iterations (default 25)
#' @param epsilon regularization parameter (default 0.1)
#' @param stop.epsilon stopping criterion, when the relative gain in likelihood is below epsilon (default 0.0001)
#' @param verbose print information (default FALSE)
#' @param no_cores number of cores (default to 1)
#' @export
#' @importFrom parallel mclapply
zinb.PCA.00 = function(datamatrix, k=2, alt.number=25, epsilon=0.1, stop.epsilon=.0001, verbose=FALSE, no_cores=1){

    n <- nrow(datamatrix)
    p <- ncol(datamatrix)

    # Initialize U and V by PCA on log(count+1) matrix
    PCA.init <- prcomp(log(1+datamatrix),center=TRUE,scale.=TRUE)
    U <- PCA.init$x[,1:k]
    V <- PCA.init$rotation[,1:k]

    # Initialize W and theta to 1
    W <- matrix(0,nrow=p,ncol=k)
    a.theta <- numeric(p)
    X.theta <- matrix(1,nrow=n) # the model is theta = exp(X.theta %*% a.theta)

    total.lik=rep(NA,alt.number)

    for (alt in 1:alt.number){
        if (verbose) {cat("Iteration ",alt,"\n",sep="")}

        # Evaluate total likelihood before alternation num alt
        total.lik[alt] <- zinb.loglik(datamatrix, exp( U %*% t(V) ), exp(X.theta %*% a.theta), U %*% t(W))
        if (verbose) {cat("log-likelihood = ",total.lik[alt],"\n",sep="")}

        # If the increase in likelihood is smaller than 0.5%, stop maximization
        if(alt>1){if(abs((total.lik[alt]-total.lik[alt-1])/total.lik[alt-1])<stop.epsilon)break}


        # Fix U, optimize in V, W and theta
        ptm <- proc.time()
        estimate <- matrix(unlist( parallel::mclapply(seq(p), function(i) {
                optim( fn=zinb.loglik.regression , gr=gradient.zinb.loglik.regression , par=c(V[i,],W[i,], a.theta[i]) , Y=datamatrix[,i] , X.mu=U , X.pi=U , X.theta=X.theta , epsilon=epsilon, control=list(fnscale=-1,trace=0) , method="BFGS")$par } , mc.cores=no_cores)) , nrow=2*k+1)
        if (verbose) {print(proc.time()-ptm)}

#        estimate <- matrix(unlist(sapply(seq(p), function(i) {
#            optim( fn=zinb.loglik.regression , gr=gradient.zinb.loglik.regression , par=c(V[i,],W[i,], a.theta[i]) , Y=datamatrix[,i] , X.mu=U , X.pi=U , X.theta=X.theta , epsilon=epsilon, control=list(fnscale=-1,trace=0) , method="BFGS")$par })) , nrow=2*k+1)

        V <- t(estimate[1:k,])
        W <- t(estimate[(k+1):(2*k),])
        a.theta <- estimate[(2*k+1),]


        if (verbose) {cat("log-likelihood = ",zinb.loglik(datamatrix, exp( U %*% t(V) ), exp(X.theta %*% a.theta), U %*% t(W)),"\n",sep="")}

        # Fix V, W, theta, optimize in U
        ptm <- proc.time()
#        estimate <- sapply(seq(n), function(i) {
#            optim( fn=zinb.loglik.regression , gr=gradient.zinb.loglik.regression , par=c(U[i,]) , Y=datamatrix[i,] , Y.mu=V , Y.pi=W , offset.theta=a.theta , epsilon=epsilon, control=list(fnscale=-1,trace=0) , method="BFGS")$par })
        estimate <- matrix(unlist( parallel::mclapply(seq(n), function(i) {
            optim( fn=zinb.loglik.regression , gr=gradient.zinb.loglik.regression , par=c(U[i,]) , Y=datamatrix[i,] , Y.mu=V , Y.pi=W , offset.theta=a.theta , epsilon=epsilon, control=list(fnscale=-1,trace=0) , method="BFGS")$par } , mc.cores=no_cores)) , nrow=k)
        U <- t(estimate)
        if (verbose) {print(proc.time()-ptm)}
    }
    zinb.result <- list(U=U,V=V,W=W,theta=exp(a.theta))
}



#' Estimation of latent factors from a count matrix
#' @param datamatrix the count data matrix (cells in rows, genes in columns)
#' @param k number of latent factors (default 2)
#' @param alt.number maximum number of iterations (default 25)
#' @param epsilon regularization parameter (default 0.1)
#' @param stop.epsilon stopping criterion, when the relative gain in likelihood is below epsilon (default 0.0001)
#' @param verbose print information (default FALSE)
#' @param no_cores number of cores (default to 1)
#' @export
#' @importFrom parallel mclapply
zinb.PCA = function(datamatrix, k=2, alt.number=25, epsilon=0.1, stop.epsilon=.0001, verbose=FALSE, no_cores=1){

    n <- nrow(datamatrix)
    p <- ncol(datamatrix)

    # Initialize U and V by SVD on log(count+1) matrix
    s <- svd(log(1+datamatrix))
    U <- s$u[,1:k] %*% sqrt(diag(s$d[1:k]))
    V <- s$v[,1:k] %*% sqrt(diag(s$d[1:k]))

    # Initialize W and theta to 1
    W <- matrix(0,nrow=p,ncol=k)
    a.theta <- numeric(p)
    X.theta <- matrix(1,nrow=n) # the model is theta = exp(X.theta %*% a.theta)

    total.lik=rep(NA,alt.number)

    for (alt in 1:alt.number){
        if (verbose) {cat("Iteration ",alt,"\n",sep="")}

        # Evaluate total likelihood before alternation num alt
        total.lik[alt] <- zinb.loglik(datamatrix, exp( U %*% t(V) ), exp(X.theta %*% a.theta), U %*% t(W)) - epsilon*(sum(U^2)+sum(V^2)+sum(W^2)+sum(a.theta^2))/2
        if (verbose) {cat("log-likelihood = ",total.lik[alt],"\n",sep="")}

        # If the increase in likelihood is smaller than 0.5%, stop maximization
        if(alt>1){if(abs((total.lik[alt]-total.lik[alt-1])/total.lik[alt-1])<stop.epsilon)break}


        # Fix U, optimize in V, W and theta
        ptm <- proc.time()
        estimate <- matrix(unlist( parallel::mclapply(seq(p), function(i) {
            optim( fn=zinb.loglik.regression , gr=gradient.zinb.loglik.regression , par=c(V[i,],W[i,], a.theta[i]) , Y=datamatrix[,i] , X.mu=U , X.pi=U , X.theta=X.theta , epsilon=epsilon, control=list(fnscale=-1,trace=0) , method="BFGS")$par } , mc.cores=no_cores)) , nrow=2*k+1)
        if (verbose) {print(proc.time()-ptm)}

        V <- t(estimate[1:k,])
        W <- t(estimate[(k+1):(2*k),])
        a.theta <- estimate[(2*k+1),]

        # Orthogonalize U, V, W
        o <- orthogonalizeJointly(U,V,W)
        U <- o$U
        V <- o$V
        W <- o$W


        if (verbose) {cat("log-likelihood = ",zinb.loglik(datamatrix, exp( U %*% t(V) ), exp(X.theta %*% a.theta), U %*% t(W))- epsilon*(sum(U^2)+sum(V^2)+sum(W^2)+sum(a.theta^2))/2,"\n",sep="")}

        # Fix V, W, theta, optimize in U
        ptm <- proc.time()

        estimate <- matrix(unlist( parallel::mclapply(seq(n), function(i) {
            optim( fn=zinb.loglik.regression , gr=gradient.zinb.loglik.regression , par=c(U[i,]) , Y=datamatrix[i,] , Y.mu=V , Y.pi=W , offset.theta=a.theta , epsilon=epsilon, control=list(fnscale=-1,trace=0) , method="BFGS")$par } , mc.cores=no_cores)) , nrow=k)
        U <- t(estimate)

        if (verbose) {print(proc.time()-ptm)}

        # Orthogonalize U, V, W
        o <- orthogonalizeJointly(U,V,W)
        U <- o$U
        V <- o$V
        W <- o$W

    }
    zinb.result <- list(U=U,V=V,W=W,theta=exp(a.theta))
}

#' Estimation of latent factors from a count matrix including estimation of the size factors
#' @param datamatrix the count data matrix (cells in rows, genes in columns)
#' @param k number of latent factors (default 2)
#' @param alt.number maximum number of iterations (default 25)
#' @param epsilon regularization parameter (default 0.1)
#' @param stop.epsilon stopping criterion, when the relative gain in likelihood is below epsilon (default 0.0001)
#' @param verbose print information (default FALSE)
#' @param no_cores number of cores (default to 1)
#' @export
#' @importFrom parallel mclapply
zinb.PCA.sf = function(datamatrix, k=2, alt.number=25, epsilon=0.1, stop.epsilon=.0001, verbose=FALSE, no_cores=1){

    n <- nrow(datamatrix)
    p <- ncol(datamatrix)

    # Initialize U and V by SVD on log(count+1) matrix
    s <- svd(log(1+datamatrix))
    U <- s$u[,1:k] %*% sqrt(diag(s$d[1:k]))
    V <- s$v[,1:k] %*% sqrt(diag(s$d[1:k]))

    #create a vector of size factors initialized as 0
    SF=rep(0,n)
    V.1=rep(1,p)

    # Initialize W and theta to 1
    W <- matrix(0,nrow=p,ncol=k)
    a.theta <- numeric(p)
    X.theta <- matrix(1,nrow=n) # the model is theta = exp(X.theta %*% a.theta)

    total.lik=rep(NA,alt.number)

    for (alt in 1:alt.number){
        if (verbose) {cat("Iteration ",alt,"\n",sep="")}

        # Evaluate total likelihood before alternation num alt
        total.lik[alt] <- zinb.loglik(datamatrix, exp( cbind(U,SF) %*% t(cbind(V,V.1)) ), exp(X.theta %*% a.theta), U %*% t(W)) - epsilon*(sum(U^2)+sum(V^2)+sum(W^2)+sum(a.theta^2))/2
        if (verbose) {cat("log-likelihood = ",total.lik[alt],"\n",sep="")}

        # If the increase in likelihood is smaller than 0.5%, stop maximization
        if(alt>1){if(abs((total.lik[alt]-total.lik[alt-1])/total.lik[alt-1])<stop.epsilon)break}


        # Fix U, optimize in V, W and theta
        ptm <- proc.time()

        estimate <- matrix(unlist( parallel::mclapply(seq(p), function(i) {
            optim( fn=zinb.loglik.regression , gr=gradient.zinb.loglik.regression , par=c(V[i,],W[i,], a.theta[i]) , Y=datamatrix[,i] , X.mu=U , offset.mu=SF, X.pi=U , X.theta=X.theta , epsilon=epsilon, control=list(fnscale=-1,trace=0) , method="BFGS")$par } , mc.cores=no_cores)) , nrow=2*k+1)

        if (verbose) {print(proc.time()-ptm)}

        V <- t(estimate[1:k,])
        W <- t(estimate[(k+1):(2*k),])
        a.theta <- estimate[(2*k+1),]
        # Orthogonalize U, V, W
        o <- orthogonalizeJointly(U,V,W)
        U <- o$U
        V <- o$V
        W <- o$W


        if (verbose) {cat("log-likelihood = ",zinb.loglik(datamatrix, exp( cbind(U,SF) %*% t(cbind(V,V.1)) ), exp(X.theta %*% a.theta), U %*% t(W))- epsilon*(sum(U^2)+sum(V^2)+sum(W^2)+sum(a.theta^2))/2,"\n",sep="")}

        # Fix V, W, theta, optimize in U
        ptm <- proc.time()

        estimate <- matrix(unlist( parallel::mclapply(seq(n), function(i) {
            optim( fn=zinb.loglik.regression , gr=gradient.zinb.loglik.regression , par=c(SF[i],U[i,]) , Y=datamatrix[i,] , X.mu=matrix(1,nrow=p), Y.mu=V , Y.pi=W , offset.theta=a.theta , epsilon=epsilon, control=list(fnscale=-1,trace=0) , method="BFGS")$par } , mc.cores=no_cores)) , nrow=k+1)

        U <- t(estimate[2:3,])
        SF <- estimate[1,]

        if (verbose) {print(proc.time()-ptm)}

        # Orthogonalize U, V, W
        o <- orthogonalizeJointly(U,V,W)
        U <- o$U
        V <- o$V
        W <- o$W

    }
    zinb.result <- list(U=U,V=V,W=W,theta=exp(a.theta),SF=SF)
}


#' Estimation of latent factors from a count matrix including estimation of the logM size factors and Pi size factors
#' @param datamatrix the count data matrix (cells in rows, genes in columns)
#' @param k number of latent factors (default 2)
#' @param alt.number maximum number of iterations (default 25)
#' @param epsilon regularization parameter (default 0.1)
#' @param stop.epsilon stopping criterion, when the relative gain in likelihood is below epsilon (default 0.0001)
#' @param verbose print information (default FALSE)
#' @param no_cores number of cores (default to 1)
#' @export
#' @importFrom parallel mclapply

zinb.PCA.sf.pif = function(datamatrix, k=2, alt.number=25, epsilon=0.1, stop.epsilon=.0001, verbose=FALSE, no_cores=1){

    n <- nrow(datamatrix)
    p <- ncol(datamatrix)

    # Initialize U and V by SVD on log(count+1) matrix
    s <- svd(log(1+datamatrix))
    U <- s$u[,1:k] %*% sqrt(diag(s$d[1:k]))
    V <- s$v[,1:k] %*% sqrt(diag(s$d[1:k]))

    #create a vector of size factors initialized as 0
    SF=rep(0,n)
    PiF=rep(0,n)
    V.1=rep(1,p)

    # Initialize W and theta to 1
    W <- matrix(0,nrow=p,ncol=k)
    a.theta <- numeric(p)
    X.theta <- matrix(1,nrow=n) # the model is theta = exp(X.theta %*% a.theta)

    total.lik=rep(NA,alt.number)

    for (alt in 1:alt.number){
        if (verbose) {cat("Iteration ",alt,"\n",sep="")}

        # Evaluate total likelihood before alternation num alt
        total.lik[alt] <- zinb.loglik(datamatrix, exp( cbind(U,SF) %*% t(cbind(V,V.1)) ), exp(X.theta %*% a.theta), cbind(U,PiF) %*% t(cbind(W,V.1))) - epsilon*(sum(U^2)+sum(V^2)+sum(W^2)+sum(a.theta^2))/2
        if (verbose) {cat("log-likelihood = ",total.lik[alt],"\n",sep="")}

        # If the increase in likelihood is smaller than 0.5%, stop maximization
        if(alt>1){if(abs((total.lik[alt]-total.lik[alt-1])/total.lik[alt-1])<stop.epsilon)break}


        # Fix U, optimize in V, W and theta
        ptm <- proc.time()

        estimate <- matrix(unlist( parallel::mclapply(seq(p), function(i) {
            optim( fn=zinb.loglik.regression , gr=gradient.zinb.loglik.regression , par=c(V[i,],W[i,], a.theta[i]) , Y=datamatrix[,i] , X.mu=U , offset.mu=SF, offset.pi=PiF, X.pi=U , X.theta=X.theta , epsilon=epsilon, control=list(fnscale=-1,trace=0) , method="BFGS")$par } , mc.cores=no_cores)) , nrow=2*k+1)

        if (verbose) {print(proc.time()-ptm)}

        V <- t(estimate[1:k,])
        W <- t(estimate[(k+1):(2*k),])
        a.theta <- estimate[(2*k+1),]

        # Orthogonalize U, V, W
        o <- orthogonalizeJointly(U,V,W)
        U <- o$U
        V <- o$V
        W <- o$W


        if (verbose) {cat("log-likelihood = ",zinb.loglik(datamatrix, exp( cbind(U,SF) %*% t(cbind(V,V.1)) ), exp(X.theta %*% a.theta), cbind(U,PiF) %*% t(cbind(W,V.1)))- epsilon*(sum(U^2)+sum(V^2)+sum(W^2)+sum(a.theta^2))/2,"\n",sep="")}

        # Fix V, W, theta, optimize in U
        ptm <- proc.time()

        estimate <- matrix(unlist( parallel::mclapply(seq(n), function(i) {
            optim( fn=zinb.loglik.regression , gr=gradient.zinb.loglik.regression , par=c(SF[i],U[i,],PiF[i]) , Y=datamatrix[i,] , X.mu=matrix(1,nrow=p), X.pi=matrix(1,nrow=p), Y.mu=V , Y.pi=W , offset.theta=a.theta , epsilon=epsilon, control=list(fnscale=-1,trace=0) , method="BFGS")$par } , mc.cores=no_cores)) , nrow=k+2)

        U <- t(estimate[2:3,])
        SF <- estimate[1,]
        PiF <-estimate[4,]

        if (verbose) {print(proc.time()-ptm)}

        # Orthogonalize U, V, W
        o <- orthogonalizeJointly(U,V,W)
        U <- o$U
        V <- o$V
        W <- o$W

    }
    zinb.result <- list(U=U,V=V,W=W,theta=exp(a.theta),SF=SF,PiF=PiF)
}

#' Estimation of latent factors from a count matrix including estimation of the logM size factors and Pi size factors
#' @param datamatrix the count data matrix (cells in rows, genes in columns)
#' @param k number of latent factors (default 2)
#' @param alt.number maximum number of iterations (default 25)
#' @param epsilon regularization parameter (default 0.1)
#' @param stop.epsilon stopping criterion, when the relative gain in likelihood is below epsilon (default 0.0001)
#' @param verbose print information (default FALSE)
#' @param no_cores number of cores (default to 1)
#' @export
#' @importFrom parallel mclapply

zinb.PCA.sf.pif.center = function(datamatrix, k=2, alt.number=25, epsilon=0.1, stop.epsilon=.0001, verbose=FALSE, no_cores=1){

    n <- nrow(datamatrix)
    p <- ncol(datamatrix)

    #create a vector of size factors initialized as 0
    SF <- rep(0,n)
    PiF <- rep(0,n)
    V.1 <- rep(1,p)
    U.1 <- rep(1,n)


    # Initialize U and V by SVD on log(count+1) matrix
    s <- svd(log(1+datamatrix))
    U <- cbind ( s$u[,1:k] %*% sqrt(diag(s$d[1:k])) )
    V <- cbind ( s$v[,1:k] %*% sqrt(diag(s$d[1:k])) )
    V.U <- matrix(apply(log(1+datamatrix),2,median),ncol=1)
    W.U <- matrix(0,nrow=p,ncol=1)

    # Initialize W and theta to 1
    W <- matrix(0,nrow=p,ncol=k)
    a.theta <- numeric(p)
    X.theta <- matrix(1,nrow=n) # the model is theta = exp(X.theta %*% a.theta)

    total.lik=rep(NA,alt.number)

    for (alt in 1:alt.number){
        if (verbose) {cat("Iteration ",alt,"\n",sep="")}

        # Evaluate total likelihood before alternation num alt
        total.lik[alt] <- zinb.loglik(datamatrix, exp( cbind(U,SF,U.1) %*% t(cbind(V,V.1,V.U)) ), exp(X.theta %*% a.theta), cbind(U,PiF,U.1) %*% t(cbind(W,V.1,W.U))) - epsilon*(sum(U^2)+sum(V^2)+sum(W^2)+sum(a.theta^2))/2
        if (verbose) {cat("log-likelihood = ",total.lik[alt],"\n",sep="")}

        # If the increase in likelihood is smaller than 0.5%, stop maximization
        if(alt>1){if(abs((total.lik[alt]-total.lik[alt-1])/total.lik[alt-1])<stop.epsilon)break}


        # Fix U, optimize in V, W and theta
        ptm <- proc.time()

        estimate <- matrix(unlist( parallel::mclapply(seq(p), function(i) {
            optim( fn=zinb.loglik.regression , gr=gradient.zinb.loglik.regression , par=c(V[i,],V.U[i],W[i,],W.U[i], a.theta[i]) , Y=datamatrix[,i] , X.mu=cbind(U,U.1) , offset.mu=SF, offset.pi=PiF, X.pi=cbind(U,U.1) , X.theta=X.theta , epsilon=epsilon, control=list(fnscale=-1,trace=0) , method="BFGS")$par } , mc.cores=no_cores)) , nrow=2*k+3)

        if (verbose) {print(proc.time()-ptm)}

        V <- t(estimate[1:k,])
        V.U[,1] <- estimate[k+1,]
        W <- t(estimate[(k+2):(2*k+1),])
        W.U <- estimate[2*k+2,]
        a.theta <- estimate[(2*k+3),]

        # Orthogonalize U, V, W
        o <- orthogonalizeJointly(U,V,W)
        U <- o$U
        V <- o$V
        W <- o$W


        if (verbose) {cat("log-likelihood = ",zinb.loglik(datamatrix, exp( cbind(U,SF,U.1) %*% t(cbind(V,V.1,V.U)) ), exp(X.theta %*% a.theta), cbind(U,PiF,U.1) %*% t(cbind(W,V.1,W.U)))- epsilon*(sum(U^2)+sum(V^2)+sum(W^2)+sum(a.theta^2))/2,"\n",sep="")}

        # Fix V, W, theta, optimize in U
        ptm <- proc.time()

        estimate <- matrix(unlist( parallel::mclapply(seq(n), function(i) {
            optim( fn=zinb.loglik.regression , gr=gradient.zinb.loglik.regression , par=c(SF[i],U[i,1:k],PiF[i]) , offset.mu=V.U, offset.pi=W.U, Y=datamatrix[i,] , X.mu=matrix(1,nrow=p), X.pi=matrix(1,nrow=p), Y.mu=V , Y.pi=W , offset.theta=a.theta , epsilon=epsilon, control=list(fnscale=-1,trace=0) , method="BFGS")$par } , mc.cores=no_cores)) , nrow=k+2)

        U <- cbind(t(estimate[2:(k+1),]))
        SF <- estimate[1,]
        PiF <-estimate[k+2,]

        if (verbose) {print(proc.time()-ptm)}

        # Orthogonalize U, V, W
        o <- orthogonalizeJointly(U,V,W)
        U <- o$U
        V <- o$V
        W <- o$W

    }
    zinb.result <- list(U=U,V=V,W=W,theta=exp(a.theta),V.U=V.U,W.U=W.U,SF=SF,PiF=PiF)
}



#' Estimation of latent factors from a count matrix
#' @param datamatrix the count data matrix (cells in rows, genes in columns)
#' @param k number of latent factors (default 2)
#' @param size.fact is logical, TRUE if size factors should be estimated
#' @param pi.fact is logical, TRUE if cell specific zero probability should be estimated (one per cell)
#' @param center is logical, TRUE if "centering is performed"
#' @param alt.number maximum number of iterations (default 25)
#' @param epsilon regularization parameter (default 0.1)
#' @param stop.epsilon stopping criterion, when the relative gain in likelihood is below epsilon (default 0.0001)
#' @param verbose print information (default FALSE)
#' @param no_cores number of cores (default to 1)
#' @export
#' @importFrom parallel mclapply
zinb.PCA.full = function(datamatrix, k=2, size.fact = TRUE, pi.fact = TRUE, center = TRUE, alt.number=25, epsilon=0.1, stop.epsilon=.0001, verbose=FALSE, no_cores=1){

    n <- nrow(datamatrix)
    p <- ncol(datamatrix)

    # initialize U and V by SVD on log(count+1) matrix
    logdata <- log1p(datamatrix)

    if (center){
        logdata <- scale(logdata, center = TRUE, scale = FALSE)
    }

    s <- svd(logdata)
    U <- s$u[,1:k] %*% sqrt(diag(s$d[1:k]))
    V <- s$v[,1:k] %*% sqrt(diag(s$d[1:k]))

    # design vectors of ones
    ones.n <- matrix(1, nrow = n, ncol = 1)
    ones.J <- matrix(1, nrow = 1, ncol = p)

    # if library size is to be estimated, initialize with total counts
    if (size.fact) {
        tc <- scale(rowSums(datamatrix), center=TRUE, scale=TRUE)
        SF <- matrix (tc, nrow = n, ncol = 1)
    } else {
        SF <- NULL
    }

    # if taking into account heterogeneity of p of zero, initialize with total number of 0 genes?
    if (pi.fact) {
        PiF <- matrix(0, nrow = n, ncol = 1)
    } else {
        PiF <- NULL
    }

    # if centering data

    if (center){
        gene.fact <- matrix(0, nrow = 1, ncol = p)

    } else {
        gene.fact <- NULL
    }

    # Initialize W and theta to 1
    W <- matrix(0,nrow=p,ncol=k)
    a.theta <- numeric(p)
    X.theta <- matrix(1,nrow=n) # the model is theta = exp(X.theta %*% a.theta)

    total.lik=rep(NA,alt.number)

    for (alt in 1:alt.number){
        if (verbose) {cat("Iteration ",alt,"\n",sep="")}

        # matrices to calculate likelihood
        logmu <- U %*% t(V)
        logitP0 <- U %*% t(W)

        if (size.fact){
            logmu <- logmu + SF %*% ones.J
        }
        if (pi.fact){
            logitP0 <- logitP0 + PiF %*% ones.J
        }
        if (center){
            logmu <- logmu + ones.n %*% gene.fact
        }

        # Evaluate total likelihood before alternation num alt
        total.lik[alt] <- zinb.loglik(datamatrix, exp( logmu ), exp(X.theta %*% a.theta), logitP0) - epsilon*(sum(U^2)+sum(V^2)+sum(W^2)+sum(a.theta^2))/2
        if (verbose) {cat("log-likelihood = ",total.lik[alt],"\n",sep="")}

        # If the increase in likelihood is smaller than 0.5%, stop maximization
        if(alt>1){
          if(abs((total.lik[alt]-total.lik[alt-1])/total.lik[alt-1])<stop.epsilon){
            break
          }
        }

        # STEP 1 : fix U, optimize in V, W and theta
        if (size.fact){
            offset.mu.step1 <- SF
        } else {
            offset.mu.step1 <- 0
        }
        if (pi.fact){
            offset.pi.step1 <- PiF
        } else {
            offset.pi.step1 <- 0
        }

        if (center == TRUE){
            X.mu.step1 <- cbind(U,ones.n)
            par0.step1 <- cbind(V,t(gene.fact),W,a.theta)
        } else {
            X.mu.step1 <- U
            par0.step1 <- cbind(V,W,a.theta)
        }


        ptm <- proc.time()
        estimate <- matrix(unlist( parallel::mclapply(seq(p), function(i) {
            optim( fn=zinb.loglik.regression , gr=gradient.zinb.loglik.regression , par=par0.step1[i,] , Y=datamatrix[,i] , X.mu=X.mu.step1 , X.pi=U , offset.mu=offset.mu.step1, offset.pi=offset.pi.step1, X.theta=X.theta , epsilon=epsilon, control=list(fnscale=-1,trace=0) , method="BFGS")$par } , mc.cores=no_cores)) , nrow=2*k+1+as.numeric(center))
        if (verbose) {print(proc.time()-ptm)}

        V <- t(estimate[1:k,])

        if (center){
            gene.fact[1,] <- estimate[k+1,]
            W <- t(estimate[(k+2):(2*k+1),])
            a.theta <- estimate[(2*k+2),]
        } else {
            W <- t(estimate[(k+1):(2*k),])
            a.theta <- estimate[(2*k+1),]
        }


        # Orthogonalize U, V, W
        o <- orthogonalizeJointly(U,V,W)
        U <- o$U
        V <- o$V
        W <- o$W


     #   if (verbose) {cat("log-likelihood = ",zinb.loglik(datamatrix, exp( U %*% t(V) ), exp(X.theta %*% a.theta), U %*% t(W))- epsilon*(sum(U^2)+sum(V^2)+sum(W^2)+sum(a.theta^2))/2,"\n",sep="")}

        # STEP 2 : fix V, W, theta, optimize in U

        par0.step2 <- U

        if (center){
            offset.mu.step2 <- t(gene.fact)
        } else {
            offset.mu.step2 <- NULL
        }

        if (size.fact){
            X.mu.step2 <- t(ones.J)
            par0.step2 <- cbind(SF,par0.step2)
        } else {
            X.mu.step2 <- NULL
        }

        if (pi.fact){
            X.pi.step2 <- t(ones.J)
            par0.step2 <- cbind(par0.step2,PiF)
        } else {
            X.pi.step2 <- NULL
        }

        ptm <- proc.time()

        estimate <- matrix(unlist( parallel::mclapply(seq(n), function(i) {
            optim( fn=zinb.loglik.regression , gr=gradient.zinb.loglik.regression , par=par0.step2[i,] , Y=datamatrix[i,] , X.mu=X.mu.step2 , X.pi=X.pi.step2 , Y.mu=V , Y.pi=W , offset.mu=offset.mu.step2 , offset.theta=a.theta , epsilon=epsilon, control=list(fnscale=-1,trace=0) , method="BFGS")$par } , mc.cores=no_cores)) , nrow=k+as.numeric(size.fact)+as.numeric(pi.fact))

        if (size.fact){
            SF <- estimate[1,]
            U <- t(estimate[2:(k+1),])
        } else {
            U <- t(estimate[1:k,])
        }

        if (pi.fact){
            PiF <- estimate[nrow(estimate),]
        }


        if (verbose) {print(proc.time()-ptm)}

        # Orthogonalize U, V, W
        o <- orthogonalizeJointly(U,V,W)
        U <- o$U
        V <- o$V
        W <- o$W

    }
    zinb.result <- list(U=U,V=V,W=W,SF=SF,PiF=PiF,center=gene.fact, theta=exp(a.theta))
}
