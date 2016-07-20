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
#' Y <- matrix(rpois(60, lambda=1), 10, 6)
#' bio <- gl(2, 3)
#' time <- rnorm(6)
#' gc <- rnorm(10)
#' m <- zinbFit(matrix(10, 10, 6), X=model.matrix(~bio + time), V=model.matrix(~gc),
#'              which_X_pi=1L, which_V_mu=1L)
#' m <- zinbInitialize(m, Y)
#' @export
#' @importFrom glmnet glmnet
#' @importFrom parallel mclapply
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
    
    P <- Y > 0
    
    if(any(rowSums(P) == 0)) {
        stop(paste0("Sample ", which(rowSums(P) == 0)[1], " has only 0 counts!"))
    }
    
    if(any(colSums(P) == 0)) {
        stop(paste0("Gene ", which(colSums(P) == 0)[1], " has only 0 counts!"))
    }
    
    Z <- 1 - P
    L <- matrix(nrow=n, ncol=J)
    L[P] <- log(Y[P]) - m@O_mu[P]
    
    for(mm in seq_len(nb.repeat)) {
        
        Xbeta_mu <- getX_mu(m) %*% m@beta_mu
        
        ## optimize gamma_mu (in parallel for each sample)
        if(NCOL(getV_mu(m)) == 0) {
            ## do nothing: no need to estimate gamma_mu
        } else if(NCOL(getV_mu(m)) == 1) {
            if(m@V_mu_intercept) {
                ## if only intercept, regular regression
                gamma_mu <- lapply(seq(n), function(i) {
                    coef(lm.fit(x = getV_mu(m)[P[i,], , drop=FALSE], 
                                y = L[i,P[i,]], 
                                offset = Xbeta_mu[i, P[i,]]))
                })
            } else {
                ## if no intercept parcor::lm.ridge.univariate adds an intercept
                ## we don't want that so I re-wrote it
                gamma_mu <- lapply(seq(n), function(i) {
                    as.vector(ridge_one(x=getV_mu(m)[P[i,], , drop=FALSE],
                                        y=L[i,P[i,]] - Xbeta_mu[i, P[i,]],
                                        epsilon = m@penalty_gamma_mu))
                })
            }
            ## I'm not sure if we need the previous if statement.
            ## In principle, when we have an intercept, its penalty will be 0.
            ## So ridge_one should accomodate all cases. But lm may be faster?
        } else {
            ## use glmnet: we don't want to penalize the intercept, so we need
            ## to treat differently the case with / without it.
            # fit <- glmnet::glmnet(x = getV_mu(m)[P[i,], , drop=FALSE],
            #                       y = L[i,P[i,]],
            #                       offset=Xbeta_mu[i, P[i,]],
            #                       intercept=m@V_mu_intercept,
            #                       standardize=FALSE,
            #                       family="gaussian")
            # gamma_mu <- lapply(seq(n), function(i) {
            #     as.vector(coef())})
        }
        m@gamma_mu <- matrix(unlist(gamma_mu), nrow=NCOL(getV_mu(m)))
        
        tVgamma_mu <- t(getV_mu(m) %*% m@gamma_mu)

        ## optimize beta_mu (in parallel for each gene)
        if(NCOL(getX_mu(m)) == 0) {
            ## do nothing: no need to estimate beta_mu
        } else if(NCOL(getX_mu(m)) == 1) {
            beta_mu <- lapply(seq(J), function(j) {
                as.vector(ridge_one(x=getX_mu(m)[P[,j], , drop=FALSE],
                                    y=L[P[,j],j] - tVgamma_mu[P[,j], j],
                                    epsilon = m@penalty_beta_mu))
            })
        } else if(NCOL(getX_mu(m, intecept=FALSE)) == 1) {
            beta_mu <- lapply(seq(J), function(j) {
                ## lm.ridge.univariate assumes intercept
                as.vector(parcor::lm.ridge.univariate(x=getX_mu(m, intercept=FALSE)[P[,j], , drop=FALSE],
                                              y=L[P[,j],j] - tVgamma_mu[P[,j], j],
                                              lambda=m@penalty_beta_mu[-1],
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
            })
        }
        m@beta_mu <- matrix(unlist(beta_mu),nrow=NCOL(getX_mu(m)))
    }
    m
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
