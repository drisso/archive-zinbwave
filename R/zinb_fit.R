#' @describeIn zinbFit Y is a matrix of counts (genes in rows).
#' @export
#' 
#' @details By default, i.e., if no arguments other than \code{Y} are passed,
#' the model is fitted with an intercept for the regression across-samples and
#' one intercept for the regression across genes, both for mu and for pi.
#' 
#' @details This means that by default the model is fitted with \code{X_mu = X_pi = 1_n}
#' and \code{V_mu = V_pi = 1_J}. If the user explicitly passes the design matrices,
#' this behavior is overwritten, i.e., the user needs to explicitly include the intercept
#' in the design matrices.
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
                   penalty_gamma_mu, penalty_gamma_pi, ...) {

    J = NROW(Y) # number of genes
    n = NCOL(Y) # number of samples

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
        val <- max(1, length(which_X_mu)*NROW(X))
        penalty_beta_mu <- rep(1/val, length(which_X_mu))
        penalty_beta_mu[which(colSums(X_mu==rep(1, n))==n)] <- 0
    }

    if(missing(penalty_beta_pi)) {
        val <- max(1, length(which_X_pi)*NROW(X))
        penalty_beta_pi <- rep(1/val, length(which_X_pi))
        penalty_beta_pi[which(colSums(X_pi==rep(1, n))==n)] <- 0
    }
    
    if(missing(penalty_gamma_mu)) {
        val <- max(1, length(which_V_mu)*NROW(V))
        penalty_gamma_mu <- rep(1/val, length(which_V_mu))
        penalty_gamma_mu[which(colSums(V_mu==rep(1, J))==J)] <- 0
    }

    if(missing(penalty_gamma_pi)) {
        val <- max(1, length(which_V_pi)*NROW(V))
        penalty_gamma_pi <- rep(1/val, length(which_V_pi))
        penalty_gamma_pi[which(colSums(V_pi==rep(1, J))==J)] <- 0
    }

    # Create a ZinbModel object
    m <- zinbModel(n=n, J=J, X=X, V=V, which_X_mu=which_X_mu,
                   which_X_pi=which_X_pi, which_V_mu=which_V_mu,
                   which_V_pi=which_V_pi, penalty_beta_mu=penalty_beta_mu,
                   penalty_beta_pi=penalty_beta_pi, penalty_gamma_mu=penalty_gamma_mu,
                   penalty_gamma_pi=penalty_gamma_pi, ...)
    
    # Initialize the parameters
    m <- zinbInitialize(m, Y)
    
    # Optimize parameters
    m <- zinbOptimize(m, Y)
    
    m
})

#' Initialize the parameters of a ZINB regression model
#'
#' The initialization performs quick optimization of the parameters with several simplifying assumptions compared to the true model: non-zero counts are models as log-Gaussian, zeros are modeled as dropouts. The dispersion parameter is not modified.
#' @param m The model of class ZinbModel
#' @param Y The matrix of counts.
#' @return An object of class ZinbModel similar to the one given as argument with modified parameters alpha_mu, alpha_pi, beta_mu, beta_pi, gamma_mu, gamma_pi, W.
#' @examples
#' Y = matrix(10,3,5)
#' m = zinbModel(n=NROW(Y),J=NCOL(Y))
#' m = zinbInitialize(m,Y)
#' @export
zinbInitialize <- function(m, Y) {

    # TODO: write this function    
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
