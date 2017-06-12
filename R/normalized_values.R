#' Log-likelihood of the zero-inflated negative binomial model for each entry
#' in the matrix of counts
#'
#' Given a matrix of counts, this function computes the
#' log-probabilities of the counts under a zero-inflated negative binomial
#' (ZINB) model. For each count, the ZINB distribution is parametrized by three
#' parameters: the mean value and the dispersion of the negative binomial
#' distribution, and the probability of the zero component.
#'
#' @param model the zinb model
#' @param x the matrix of counts
#' @return the matrix of log-likelihood of the model.
#' @importFrom stats dnbinom
zinb.loglik.matrix <- function(model, x) {
    mu <- t(getMu(model))
    theta <- getTheta(model)
    theta_mat <- matrix(rep(theta, ncol(x), ncol = ncol(x)))
    pi <- t(getPi(model))
    log( pi * (x == 0) + (1 - pi) * dnbinom(x, size = theta, mu = mu) )
}


#' Deviance residuals of the zero-inflated negative binomial model
#'
#' Given a matrix of counts, this function computes the
#' deviance residuals under a zero-inflated negative binomial
#' (ZINB) model.
#'
#' @param model the zinb model
#' @param x the matrix of counts
#' @param ignoreW logical, if true matrix \code{W} is ignored. Default is TRUE.
#' @export
#' @return the matrix of deviance residuals of the model.
computeDevianceResiduals <- function(model, x, ignoreW = TRUE) {
    if (ignoreW) model@W <- matrix(0, ncol = nFactors(model), nrow = nSamples(model))
    mu_hat <- t(getMu(model))
    pi_hat <- t(getPi(model))
    x_hat <- (1 - pi_hat) * mu_hat
    ll <- zinb.loglik.matrix(model, x)
    sign <- 1*(x - x_hat > 0)
    sign[sign == 0] <- -1
    sign * sqrt(-2 * ll)
}



#' @describeIn zinbDimRed Y is a
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#' @export
#'
#' @param X The design matrix containing sample-level covariates, one sample per
#'   row. If missing, X will contain only an intercept. If Y is a
#'   SummarizedExperiment object, X can be a formula using the variables in the
#'   colData slot of Y.
#' @param V The design matrix containing gene-level covariates, one gene
#'   per row. If missing, V will contain only an intercept. If Y is a
#'   SummarizedExperiment object, V can be a formula using the variables in the
#'   rowData slot of Y.
#' @param commondispersion Whether or not a single dispersion for all features
#'   is estimated (default TRUE).
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
#' @param verbose Print helpful messages.
#' @param nb.repeat.initialize Number of iterations for the initialization of
#'   beta_mu and gamma_mu.
#' @param maxiter.optimize maximum number of iterations for the optimization
#'   step (default 25).
#' @param stop.epsilon.optimize stopping criterion in the optimization step,
#'   when the relative gain in likelihood is below epsilon (default 0.0001).
#' @param normalizedValues indicates wether or not you want to compute
#' normalized values for the counts after adjusting for gene and cell-level
#' covariates.
#' @param residuals indicates wether or not you want to compute the residuals
#' of the ZINB model. Deviance residuals are computed.
#'
#' @details For visualization (heatmaps, ...), please use the normalized values.
#' It corresponds to the deviance residuals when the \code{W} is not included
#' in the model but the gene and cell-level covariates are. As a results, when
#' \code{W} is not included in the model, the deviance residuals should capture
#' the biology.
#'
#' @import SummarizedExperiment
#'
#' @examples
#' se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
#'                            colData = data.frame(bio = gl(2, 3)))
#'
#' m <- zinbDimRed(se, X=model.matrix(~bio, data=colData(se)))
setMethod("zinbDimRed", "SummarizedExperiment",
          function(Y, X, V, commondispersion=TRUE, verbose=FALSE,
                   nb.repeat.initialize=2, maxiter.optimize=25,
                   stop.epsilon.optimize=.0001,
                   BPPARAM=BiocParallel::bpparam(),
                   normalizedValues = TRUE, residuals = FALSE, ...) {

              res <- zinbFit(Y, X, V, commondispersion,
                             verbose, nb.repeat.initialize, maxiter.optimize,
                             stop.epsilon.optimize, BPPARAM, ...)


              # Returns a summarizedExperiment object where normalized values
              # and deviance residuals can be added to the list of assays and
              # the W has been added to the colData matrix if K > 0
              if (normalizedValues){
                  norm <- computeDevianceResiduals(res, assay(Y), ignoreW = T)
                  assays(Y)[['normalizedValues']] <- norm
              }

              if (residuals){
                  devres <- computeDevianceResiduals(res, assay(Y), ignoreW = F)
                  assays(Y)[['residuals']] <- devres
              }

              if (nFactors(res) > 0){
                  W <- data.frame(getW(res))
                  colnames(W) <- paste0('W', seq_len(nFactors(res)))
                  colData(Y) <- cbind(colData(Y), W)
              }

              return(Y)
          }
)







