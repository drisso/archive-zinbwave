#' @describeIn zinbFit Y is a matrix of counts (genes in rows).
#' @export
#' 
#' @details By default, i.e., if no arguments other than \code{Y} are passed, 
#'   the model is fitted with an intercept for the regression across-samples and
#'   one intercept for the regression across genes, both for mu and for pi.
#'   
#' @details This means that by default the model is fitted with \code{X_mu = 
#'   X_pi = 1_n} and \code{V_mu = V_pi = 1_J}. If the user explicitly passes the
#'   design matrices, this behavior is overwritten, i.e., the user needs to 
#'   explicitly include the intercept in the design matrices.
#'   
#' @seealso \code{\link[stats]{model.matrix}}.
#' 
#' @examples 
#' bio <- gl(2, 3)
#' m <- zinbFit(matrix(10, 10, 6), X=model.matrix(~bio))
setMethod("zinbFit", "matrix",
          function(Y, X=matrix(1, ncol=1, nrow=NCOL(Y)), 
                   V=matrix(1, ncol=1, nrow=NROW(Y)), 
                   which_X_mu=seq_len(NCOL(X)), which_V_mu=seq_len(NCOL(V)),
                   which_X_pi=seq_len(NCOL(X)), which_V_pi=seq_len(NCOL(V)),
                   penalty_beta_mu, penalty_beta_pi,
                   penalty_gamma_mu, penalty_gamma_pi, ncores=1, ...) {

    J <- NROW(Y) # number of genes
    n <- NCOL(Y) # number of samples

    ## Check singularity of X and V
    X_mu <- X[, which_X_mu, drop=FALSE]
    X_pi <- X[, which_X_pi, drop=FALSE]
    V_mu <- V[, which_V_mu, drop=FALSE]
    V_pi <- V[, which_V_pi, drop=FALSE]

    if(qr(X_mu)$rank < NCOL(X_mu)) {
       stop("The matrix X_mu is singular!")
    }
    if(qr(X_pi)$rank < NCOL(X_pi)) {
        stop("The matrix X_pi is singular!")
    }
    if(qr(V_mu)$rank < NCOL(V_mu)) {
        stop("The matrix V_mu is singular!")
    }
    if(qr(V_pi)$rank < NCOL(V_pi)) {
        stop("The matrix V_pi is singular!")
    }

    ## Identify intercept and put its penalty to 0
    if(missing(penalty_beta_mu)) {
        idx <- which(colSums(X_mu == rep(1, n)) == n)
        val <- max(1, (length(which_X_mu) - length(idx)) * J)
        penalty_beta_mu <- rep(1/val, length(which_X_mu))
        penalty_beta_mu[idx] <- 0
    }

    if(missing(penalty_beta_pi)) {
        idx <- which(colSums(X_pi == rep(1, n)) == n)
        val <- max(1, (length(which_X_pi) - length(idx)) * J)
        penalty_beta_pi <- rep(1/val, length(which_X_pi))
        penalty_beta_pi[idx] <- 0
    }
    
    if(missing(penalty_gamma_mu)) {
        idx <- which(colSums(V_mu == rep(1, J)) == J)
        val <- max(1, (length(which_V_mu) - length(idx)) * n)
        penalty_gamma_mu <- rep(1/val, length(which_V_mu))
        penalty_gamma_mu[idx] <- 0
    }

    if(missing(penalty_gamma_pi)) {
        idx <- which(colSums(V_pi == rep(1, J)) == J)
        val <- max(1, (length(which_V_pi) - length(idx)) * n)
        penalty_gamma_pi <- rep(1/val, length(which_V_pi))
        penalty_gamma_pi[idx] <- 0
    }

    # Create a ZinbModel object
    m <- zinbModel(n=n, J=J, X=X, V=V, which_X_mu=which_X_mu,
                   which_X_pi=which_X_pi, which_V_mu=which_V_mu,
                   which_V_pi=which_V_pi, penalty_beta_mu=penalty_beta_mu,
                   penalty_beta_pi=penalty_beta_pi, penalty_gamma_mu=penalty_gamma_mu,
                   penalty_gamma_pi=penalty_gamma_pi, ...)
    
    # Initialize the parameters
    m <- zinbInitialize(m, Y, ncores = ncores)
    
    # Optimize parameters
    m <- zinbOptimize(m, Y)
    
    validObject(m)
    m
})

#' Initialize the parameters of a ZINB regression model
#' 
#' The initialization performs quick optimization of the parameters with several
#' simplifying assumptions compared to the true model: non-zero counts are 
#' models as log-Gaussian, zeros are modeled as dropouts. The dispersion 
#' parameter is not modified.
#' @param m The model of class ZinbModel
#' @param Y The matrix of counts.
#' @param nb.repeat Number of iterations for the estimation of beta_mu and
#'   gamma_mu.
#' @param ncores The number of cores. To be passed to mclapply.
#' @return An object of class ZinbModel similar to the one given as argument 
#'   with modified parameters alpha_mu, alpha_pi, beta_mu, beta_pi, gamma_mu, 
#'   gamma_pi, W.
#' @examples
#' Y <- matrix(rpois(60, lambda=2), 10, 6)
#' bio <- gl(2, 3)
#' time <- rnorm(6)
#' gc <- rnorm(10)
#' m <- zinbFit(Y, X=model.matrix(~bio + time), V=model.matrix(~gc),
#'              which_X_pi=1L, which_V_mu=1L, K=1)
#' m <- zinbInitialize(m, Y)
#' @export
#' @importFrom glmnet glmnet
#' @importFrom parallel mclapply
#' @importFrom softImpute softImpute
zinbInitialize <- function(m, Y, nb.repeat=2, ncores=1) {

    ## we want to work with genes in columns
    Y <- t(Y)
    
    n <- NROW(Y)
    J <- NCOL(Y)
    
    if(n != nSamples(m)) {
        stop(paste0("Y needs to have ", nSamples(m), " rows (genes)"))
    }

    if(J != nFeatures(m)) {
        stop(paste0("Y needs to have ", nFeatures(m), " columns (samples)"))
    }
    
    ## 1. Define P
    P <- Y > 0
    
    if(any(rowSums(P) == 0)) {
        stop(paste0("Sample ", which(rowSums(P) == 0)[1], " has only 0 counts!"))
    }
    
    if(any(colSums(P) == 0)) {
        stop(paste0("Gene ", which(colSums(P) == 0)[1], " has only 0 counts!"))
    }
    
    ## 2. Define L
    L <- matrix(nrow=n, ncol=J)
    L[P] <- log(Y[P]) - m@O_mu[P]

    ## 3. Define Z
    Z <- 1 - P
    
    ## 4. Estimate gamma_mu and beta_mu
    for(iter in seq_len(nb.repeat)) {
        
        Xbeta_mu <- getX_mu(m) %*% m@beta_mu
        
        ## optimize gamma_mu (in parallel for each sample)
        if(NCOL(getV_mu(m)) == 0) {
            ## do nothing: no need to estimate gamma_mu
            gamma_mu <- as.vector(m@gamma_mu)
        } else if(NCOL(getV_mu(m)) == 1) {
            ## if no intercept parcor::lm.ridge.univariate adds an intercept
            ## we don't want that so I re-wrote it
            gamma_mu <- lapply(seq(n), function(i) {
                as.vector(ridge_one(x=getV_mu(m)[P[i,], , drop=FALSE],
                                    y=L[i,P[i,]] - Xbeta_mu[i, P[i,]],
                                    epsilon = getEpsilon_gamma_mu(m)))
            })
        } else if(NCOL(getV_mu(m, intercept=FALSE)) == 1) {
            ## lm.ridge.univariate assumes intercept
            gamma_mu <- lapply(seq(n), function(i) {
                as.vector(parcor::lm.ridge.univariate(x=getV_mu(m, intercept=FALSE)[P[i,], , drop=FALSE],
                                                      y=L[i,P[i,]] - Xbeta_mu[i, P[i,]],
                                                      lambda=getEpsilon_gamma_mu(m)[-1],
                                                      scale=FALSE))
            })
        } else {
            ## use glmnet
            gamma_mu <- lapply(seq(n), function(i) {
                fit <- glmnet::glmnet(x = getV_mu(m, intercept=FALSE)[P[i,], , drop=FALSE],
                                      y = L[i,P[i,]],
                                      alpha=0,
                                      offset=Xbeta_mu[i, P[i,]],
                                      intercept=m@V_mu_intercept,
                                      standardize=FALSE,
                                      family="gaussian")
                
                ## for now, assuming all penalty_gamma_mu are the same
                as.vector(coef(fit, s=getEpsilon_gamma_mu(m)[length(m@penalty_gamma_mu)]))
            })
        }
        m@gamma_mu <- matrix(unlist(gamma_mu), nrow=NCOL(getV_mu(m)))
        
        tVgamma_mu <- t(getV_mu(m) %*% m@gamma_mu)

        ## optimize beta_mu (in parallel for each gene)
        if(NCOL(getX_mu(m)) == 0) {
            ## do nothing: no need to estimate beta_mu
            beta_mu <- as.vector(m@beta_mu)
        } else if(NCOL(getX_mu(m)) == 1) {
            beta_mu <- lapply(seq(J), function(j) {
                as.vector(ridge_one(x=getX_mu(m)[P[,j], , drop=FALSE],
                                    y=L[P[,j],j] - tVgamma_mu[P[,j], j],
                                    epsilon = getEpsilon_beta_mu(m)))
            })
        } else if(NCOL(getX_mu(m, intercept=FALSE)) == 1) {
            beta_mu <- lapply(seq(J), function(j) {
                as.vector(parcor::lm.ridge.univariate(x=getX_mu(m, intercept=FALSE)[P[,j], , drop=FALSE],
                                              y=L[P[,j],j] - tVgamma_mu[P[,j], j],
                                              lambda=getEpsilon_beta_mu(m)[-1],
                                              scale=FALSE))
            })
        } else {
            beta_mu <- lapply(seq(J), function(j) {
                fit <- glmnet::glmnet(x = getX_mu(m, intercept=FALSE)[P[,j], , drop=FALSE],
                                      y = L[P[,j],j],
                                      alpha=0,
                                      offset=tVgamma_mu[P[,j], j],
                                      intercept=m@X_mu_intercept,
                                      standardize=FALSE,
                                      family="gaussian")

                ## for now, assuming all penalty_beta_mu are the same
                as.vector(coef(fit, s=getEpsilon_beta_mu(m)[length(m@penalty_beta_mu)]))
            })
        }
        m@beta_mu <- matrix(unlist(beta_mu), nrow=NCOL(getX_mu(m)))
    }
    
    ## 5. Estimate W and alpha (only if K>0)
    if(nFactors(m) > 0) {
        ## Define D
        D <- L - getX_mu(m) %*% m@beta_mu - t(getV_mu(m) %*% m@gamma_mu)
        ## to follow JP's document, we should use type=svd, but for speedup the
        ## man page suggests the default als
        R <- softImpute::softImpute(D, 
                                    lambda=sqrt(getEpsilon_W(m) * getEpsilon_alpha_mu(m)),
                                    rank.max=nFactors(m))
        m@W <- (getEpsilon_alpha_mu(m)/getEpsilon_W(m))^(1/4) * R$u %*% diag(sqrt(R$d), nrow = length(R$d))
        m@alpha_mu <- (getEpsilon_W(m)/getEpsilon_alpha_mu(m))^(1/4) * diag(sqrt(R$d), nrow = length(R$d)) %*% t(R$v)
    }
    
    ## 6. Estimate beta_pi, gamma_pi, and alpha_pi
    XW <- cbind(getX_pi(m), m@W)
    betaalpha <- matrix(0, nrow=ncol(XW), ncol=J)
    XWbetaalpha <- XW %*% betaalpha
    XW_intercept <- all(XW[,1] == 1)
    XW_noint <- if(XW_intercept) XW[,-1, drop=FALSE] else XW
    
    for (iter in seq_len(nb.repeat)) {
        ## optimize gamma_pi
        if(NCOL(getV_pi(m)) == 0) {
            ## do nothing: no need to estimate gamma_pi
            gamma_pi <- as.vector(m@gamma_pi)
        } else if(NCOL(getV_pi(m)) == 1) {
            if(m@V_pi_intercept) {
                ## simply set gamma_pi to the logit of the proportion of zeros in each cell
                gamma_pi <- qlogis(rowSums(Z)/NCOL(Z))
            } else {
                ## TODO: implement univariate w/o intercept case
                gamma_pi <- as.vector(m@gamma_pi)
            }
        } else if(NCOL(getV_pi(m, intercept=FALSE)) == 1) {
            ## TODO: implement univariate w/ intercept case
            gamma_pi <- as.vector(m@gamma_pi)
        } else {
            gamma_pi <- lapply(seq(n), function(i) {
                fit <- glmnet::glmnet(getV_pi(m, intercept=FALSE),
                                      Z[i,],
                                      offset = XWbetaalpha[i,],
                                      alpha=0,
                                      family="binomial",
                                      standardize=FALSE,
                                      intercept=m@V_pi_intercept)
                ## for now, assuming all penalty_gamma_pi are the same
                as.vector(coef(fit, s=getEpsilon_gamma_pi(m)[length(m@penalty_gamma_pi)]))
            })
        }
        
        m@gamma_pi <- matrix(unlist(gamma_pi), nrow=NCOL(getV_pi(m)))
        
        tVgamma_pi <- t(getV_pi(m) %*% m@gamma_pi)
        
        ## optimize alpha_pi and beta_pi
        if(NCOL(XW) == 0) {
            ## do nothing: no need to estimate beta_pi and alpha_pi
        } else if(NCOL(XW) == 1) {
            if(XW_intercept) {
                ## TODO: implement intercept only case
            } else {
                ## TODO: implement univariate w/o intercept case
            }
        } else if(NCOL(XW) == 2 && XW_intercept) {
            ## TODO: implement univariate w/ intercept case
        } else {
            betaalpha_pi <- lapply(seq(J), function(j) {
                fit <- glmnet::glmnet(XW_noint,
                                      Z[,j],
                                      offset = tVgamma_pi[,j],
                                      alpha=0,
                                      family="binomial",
                                      standardize=FALSE,
                                      intercept=XW_intercept)
                
                ## for now, assuming all penalty_beta_pi are the same
                beta_pi <- as.vector(coef(fit, s=getEpsilon_beta_pi(m)[length(m@penalty_beta_pi)]))[seq_len(NCOL(m@X))]
                alpha_pi <- as.vector(coef(fit, s=getEpsilon_alpha_pi(m)[length(m@penalty_alpha_pi)]))[(NCOL(m@X)+1):NCOL(XW)]
                c(beta_pi, alpha_pi)
            })
            m@beta_pi <- matrix(unlist(lapply(betaalpha_pi, function(x) x[seq_len(NCOL(m@X))])), nrow=NCOL(getX_pi(m)))
            m@alpha_pi <- matrix(unlist(lapply(betaalpha_pi, function(x) x[(NCOL(m@X)+1):NCOL(XW)])), nrow=NCOL(m@W))
        }
    }
    
    ## 7. Initialize dispersion to 1
    m@logtheta <- rep(0, J)
    
    validObject(m)
    
    return(m)
}

#' Optimize the parameters of a ZINB regression model
#'
#' The parameters of the model given as argument are optimized by penalized maximum likelihood on the count matrix given as argument. It is recommended to call zinb_initialize before this function to have good starting point for optimization, since the optimization problem is not convex and can only converge to a local minimum.
#' @param m The model of class ZinbModel
#' @param Y The matrix of counts.
#' @return An object of class ZinbModel similar to the one given as argument with modified parameters alpha_mu, alpha_pi, beta_mu, beta_pi, gamma_mu, gamma_pi, W.
#' @examples
#' Y = matrix(10,3,5)
#' m = zinbModel(n=NROW(Y),J=NCOL(Y))
#' m = zinbInitialize(m,Y)
#' m = zinbOptimize(m,Y)
#' @export
zinbOptimize <- function(m, Y) {
    
    # TODO: write this function    
    m
}

## helper functions
ridge_one <- function(x, y, epsilon=rep(0, NCOL(x))) {
    solve(crossprod(x) + tcrossprod(epsilon)) %*% crossprod(x, y)
}
