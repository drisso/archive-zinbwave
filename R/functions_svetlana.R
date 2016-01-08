
# TODO: dangerous to define global variables since they may be overriden by other parts. Better to put them as function arguments when you need them

linkstr <- "logit"
linkobj <- make.link(linkstr)
linkinv <- linkobj$linkinv


# TODO: to document a function, do as in the 'standard.R' file of Davide, so that the documentation will be automatically created by roxygen 


#' Log-likelihood function of the zero-inflated negative binomial model for the estimation of "right parts" with 
#' log(M)=X.M*alpha.M+U*V; logit(Pi)=X.pi*alpha.pi+U*W; M=E[Y] and Y data matrix
#' @param parms a vector of parameters: j-th column of alpha.M, followed by j-th column of V, followed by j-th column of alpha.pi, followed by j-th column of W, followed by the log(1/phi)
#' @param j is the number of column (gene) to be optimized over
#' @param Y the data matrix (cells in rows, genes in columns)
#' @param X.M is the design matrix for M 
#' @param X.pi is the design matrix for Pi.
#' @param U is the current value for the latent factor matrix
#' @param offsetx the offset for the regression on M
#' @param offsetw the offset for the regression on Pi
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())


ziNegBin <- function(parms,j,X.M,X.pi,U,Y,epsilon=0,offsetx=0,offsetz=0) {
    X=cbind(X.M,U)
    Z=cbind(X.pi,U)
    kx=ncol(X)
    kz=ncol(Z)
    # TODO : I added 'epsilon=0' as parameter otherwise it was undefined. Check this is correct
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    phi <- as.vector(linkinv(Z %*% parms[(kx + 1):(kx + kz)] + 
                                 offsetz))
    theta <- exp(parms[(kx + kz) + 1])
    loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(0, 
                                                                     size = theta, mu = mu, log = TRUE))))
    loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(Y[,j], 
                                                       size = theta, mu = mu, log = TRUE))
    Y0 <- Y <= 0 #==1 if counts is 0
    Y1 <- Y > 0 #==1 if count is not 0
    loglik <- sum(loglik0[Y0[,j]]) + sum(loglik1[Y1[,j]])
    loglik-sum(parms[(kx + 1):(kx + kz)]^2)*epsilon
}

#' Gradient with respect to "right parts" of log-likelihood function of the zero-inflated negative binomial model with 
#' log(M)=X.M*alpha.M+U*V; logit(Pi)=X.pi*alpha.pi+U*W; M=E[Y] and Y data matrix
#' @param parms a vector of parameters: j-th column of alpha.M, followed by j-th column of V, followed by j-th column of alpha.pi, followed by j-th column of W, followed by the log(1/phi)
#' @param j is the number of column (gene) to be optimized over
#' @param Y the data matrix (cells in rows, genes in columns)
#' @param X.M is the design matrix for M 
#' @param X.pi is the design matrix for Pi.
#' @param U is the current value for the latent factor matrix
#' @param offsetx the offset for the regression on M
#' @param offsetw the offset for the regression on Pi
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())

#gradient function
gradNegBin <- function(parms,j,X.M,X.pi,U,Y,epsilon=0,offsetx=0,offsetz=0) {
    X=cbind(X.M,U)
    Z=cbind(X.pi,U)
    kx=ncol(X)
    kz=ncol(Z)
    eta <- as.vector(X %*% parms[1:kx] + offsetx)
    mu <- exp(eta)
    etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
    muz <- linkinv(etaz)
    theta <- exp(parms[(kx + kz) + 1])
    Y0 <- Y <= 0 #==1 if counts is 0
    Y1 <- Y > 0 #==1 if count is not 0
    clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
    dens0 <- muz * (1 - as.numeric(Y1[,j])) + exp(log(1 - muz) + 
                                                      clogdens0)
    wres_count <- ifelse(Y1[,j], Y[,j] - mu * (Y[,j] + theta)/(mu + theta), 
                         -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(theta) - 
                                  log(mu + theta) + log(mu)))
    wres_zero <- ifelse(Y1[,j], -1/(1 - muz) * linkobj$mu.eta(etaz), 
                        (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
    wres_theta <- theta * ifelse(Y1[,j], digamma(Y[,j] + theta) - 
                                     digamma(theta) + log(theta) - log(mu + theta) + 1 - 
                                     (Y[,j] + theta)/(mu + theta), exp(-log(dens0) + log(1 - 
                                                                                             muz) + clogdens0) * (log(theta) - log(mu + theta) + 
                                                                                                                      1 - theta/(mu + theta)))
    colSums(cbind(wres_count * X, wres_zero * Z-2*parms[(kx + 1):(kx + kz)]*epsilon, wres_theta))
}



######################################################################################
#Optimization with respect to U
######################################################################################

#' Log-likelihood function of the zero-inflated negative binomial model for the estimation of U in the model 
#' log(M)=X.M*alpha.M+U*V; logit(Pi)=X.pi*alpha.pi+U*W; M=E[Y] with Y data matrix
#' @param parms a vector of parameters: i-th line of U
#' @param i is the number of ligne (cell) to be optimized over
#' @param Y the data matrix (cells in rows, genes in columns)
#' @param X.M is the design matrix for M 
#' @param X.pi is the design matrix for Pi.
#' @param V,W,alpha.M,alpha.pi are the "right parts" to be given as arguments
#' @theta vector of length J with gene specific dispersion parameters
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())

#log(M)=X.M%*%alpa.M+U%*%V, transpose the equation and consider X.M%*%alpha.M as known

ziNegBin.U <- function(parms,i,V,W,alpha.M,alpha.pi,Y,theta){
    offset.M=t(alpha.M)%*%t(X.M)
    offset.pi=t(alpha.pi)%*%t(X.pi)
    theta=exp(theta)
    Y0 <- Y == 0
    Y1 <- Y != 0
    mu <- as.vector(exp(t(V) %*% parms + offset.M[,i]))
    phi <- as.vector(linkinv(t(W) %*% parms + 
                                 offset.pi[,i]))
    loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(0, 
                                                                     size = theta, mu = mu, log = TRUE))))
    # counts[i,] is i-th column of t(counts)
    loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(counts[i,], 
                                                       size = theta, mu = mu, log = TRUE))
    loglik <- sum(loglik0[Y0[i,]]) + sum(loglik1[Y1[i,]])
    loglik
}
#sum(sapply(seq(1000),function(i){ziNegBin.U(U[i,],i,V,W,t(alpha.M)%*%t(X.M),t(alpha.pi)%*%t(X.pi),counts,thetas)}))
# gradient function for optimization in U
#sum(sapply(seq(20),function(j){ziNegBin(trueparam[,j],j)}))
# gradient function for optimization in U

#' Gradient of the log-likelihood function of the zero-inflated negative binomial model for the estimation of U in the model 
#' log(M)=X.M*alpha.M+U*V; logit(Pi)=X.pi*alpha.pi+U*W; M=E[Y] with Y data matrix
#' @param parms a vector of parameters: i-th line of U
#' @param i is the number of ligne (cell) to be optimized over
#' @param Y the data matrix (cells in rows, genes in columns)
#' @param X.M is the design matrix for M 
#' @param X.pi is the design matrix for Pi.
#' @param V,W,alpha.M,alpha.pi are the "right parts" to be given as arguments
#' @theta vector of length J with gene specific dispersion parameters
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())

gradNegBin.U <- function(parms,i,V,W,alpha.M,alpha.pi,Y,theta) {
    Y0 <- Y == 0
    Y1 <- Y != 0
    offset.M=t(alpha.M)%*%t(X.M)
    offset.pi=t(alpha.pi)%*%t(X.pi)
    eta <- as.vector(t(V) %*% parms + offset.M[,i])
    mu <- exp(eta)
    etaz <- as.vector(t(W) %*% parms + offset.pi[,i])
    muz <- linkinv(etaz)
    theta <- exp(theta)
    clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
    dens0 <- muz * (1 - as.numeric(Y1[i,])) + exp(log(1 - muz) + 
                                                            clogdens0)
    wres_count <- ifelse(Y1[i,],  counts[i,] - mu * (counts[i,] + theta)/(mu + thetas), 
                         -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(theta) - 
                                  log(mu + theta) + log(mu)))
    wres_zero <- ifelse(Y1[i,], -1/(1 - muz) * linkobj$mu.eta(etaz), 
                        (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
    colSums(wres_count * t(V) + wres_zero * t(W))
}

