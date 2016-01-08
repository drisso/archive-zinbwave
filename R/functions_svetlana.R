
# TODO: dangerous to define global variables since they may be overriden by other parts. Better to put them as function arguments when you need them

linkstr <- "logit"
linkobj <- make.link(linkstr)
linkinv <- linkobj$linkinv


# TODO: to document a function, do as in the 'standard.R' file of Davide, so that the documentation will be automatically created by roxygen 



#loglikelihood
#parms should be a vector of parameters for gene j
#first all M-parameters, then all Pi-parameters, then theta

ziNegBin <- function(parms,j,X,Z,Y,epsilon=0) {
    # TODO : I added 'epsilon=0' as parameter otherwise it was undefined. Check this is correct
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    # TODO (ERROR) : kx and kz undefined. Either an argument of the function, or define it in the function
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

#gradient function
gradNegBin <- function(parms,j,X,Z,Y,epsilon=0) {
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
    
    # TODO: I think the use of epsilon below is wrong, (epsilon!=0) is binary, why not use epsilon directly?
    colSums(cbind(wres_count * X, wres_zero * Z-2*parms[(kx + 1):(kx + kz)]*(epsilon!=0), wres_theta))
}



######################################################################################
#Optimization with respect to U
#log(M)=X.M%*%alpa.M+U%*%V, transpose the equation and consider X.M%*%alpha.M as known


# parms here are elements of i-th column of t(U) (i-th line of U) 
# all the other arguments are known
# i is the number of line of U to be optimized over
# V,W are from previous iteration
# offset.M = t(alpha.M)%*%t(X.M)
# offset.pi = t(alpha.pi)%*%t(X.pi)
# counts is the matrix of counts
# thetas is vector of thetas of all genes 

#loglik which will take U as parameter

ziNegBin.U <- function(parms,i,V,W,offset.M,offset.pi,counts,thetas){
    thetas=exp(thetas)
    is0counts <- counts == 0
    n0counts <- counts != 0
    mu <- as.vector(exp(t(V) %*% parms + offset.M[,i]))
    phi <- as.vector(linkinv(t(W) %*% parms + 
                                 offset.pi[,i]))
    loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(0, 
                                                                     size = thetas, mu = mu, log = TRUE))))
    # counts[i,] is i-th column of t(counts)
    loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(counts[i,], 
                                                       size = thetas, mu = mu, log = TRUE))
    loglik <- sum(loglik0[is0counts[i,]]) + sum(loglik1[n0counts[i,]])
    loglik
}
#sum(sapply(seq(1000),function(i){ziNegBin.U(U[i,],i,V,W,t(alpha.M)%*%t(X.M),t(alpha.pi)%*%t(X.pi),counts,thetas)}))
# gradient function for optimization in U
#sum(sapply(seq(20),function(j){ziNegBin(trueparam[,j],j)}))
# gradient function for optimization in U

gradNegBin.U <- function(parms,i,V,W,offset.M,offset.pi,counts,thetas) {
    is0counts <- counts == 0
    n0counts <- counts != 0
    eta <- as.vector(t(V) %*% parms + offset.M[,i])
    mu <- exp(eta)
    etaz <- as.vector(t(W) %*% parms + offset.pi[,i])
    muz <- linkinv(etaz)
    thetas <- exp(thetas)
    clogdens0 <- dnbinom(0, size = thetas, mu = mu, log = TRUE)
    dens0 <- muz * (1 - as.numeric(n0counts[i,])) + exp(log(1 - muz) + 
                                                            clogdens0)
    wres_count <- ifelse(n0counts[i,],  counts[i,] - mu * (counts[i,] + thetas)/(mu + thetas), 
                         -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(thetas) - 
                                  log(mu + thetas) + log(mu)))
    wres_zero <- ifelse(n0counts[i,], -1/(1 - muz) * linkobj$mu.eta(etaz), 
                        (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
    colSums(wres_count * t(V) + wres_zero * t(W))
}

