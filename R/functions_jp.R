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
    lognorm <- -logitP0 - copula::log1pexp(-logitP0)

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
            logitPi <- offset.mu
        } else {
            logitPi <- logitPi + offset.mu
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


