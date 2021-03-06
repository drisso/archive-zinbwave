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
    mu <- getMu(model)
    theta <- getTheta(model)
    theta_mat <- matrix(rep(theta, each = nrow(x)), ncol = ncol(x))
    pi <- getPi(model)
    lik <- pi * (x == 0) + (1 - pi) * dnbinom(x, size = theta_mat, mu = mu)
    lik[lik == 0] <- min(lik[lik != 0]) #to avoid log lik to be infinite
    log(lik)
}


#' Deviance residuals of the zero-inflated negative binomial model
#'
#' Given a matrix of counts, this function computes the
#' deviance residuals under a zero-inflated negative binomial
#' (ZINB) model.
#'
#' @param model the zinb model
#' @param x the matrix of counts n cells by J genes
#' @param ignoreW logical, if true matrix \code{W} is ignored. Default is TRUE.
#' @export
#' @return the matrix of deviance residuals of the model.
computeDevianceResiduals <- function(model, x, ignoreW = TRUE) {

    # this makes a copy of "model" -- is there a more efficient way?
    if (ignoreW) {
        model@W <- matrix(0, ncol = nFactors(model), nrow = nSamples(model))
    }

    mu_hat <- getMu(model)
    pi_hat <- getPi(model)
    x_hat <- (1 - pi_hat) * mu_hat
    ll <- zinb.loglik.matrix(model, x)
    sign <- 1*(x - x_hat > 0)
    sign[sign == 0] <- -1
    sign * sqrt(-2 * ll)
}

#' Impute the zeros using the estimated parameters from the ZINB model.
#'
#' Given a matrix of counts and a zinb model, this function computes the
#' imputed counts under a zero-inflated negative binomial
#' (ZINB) model.
#'
#' @param model the zinb model
#' @param x the matrix of counts n cells by J genes
#' @export
#' @return the matrix of imputed counts.
imputeZeros <- function(model, x) {
    mu <- getMu(model)
    pi <- getPi(model)
    phi <- getPhi(model)
    imputedZero <- (mu * pi) / (pi + (1 - pi) * t((1 + t(mu) * phi)^(- 1/ phi)))
    x[x == 0] <- imputedZero[x == 0]
    x
}



#' @describeIn zinbwave Y is a
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
#' @param fitted_model a \code{\link{ZinbModel}} object.
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
#' @param imputedValues indicates wether or not you want to compute the imputed
#' counts of the ZINB model.
#'
#' @details For visualization (heatmaps, ...), please use the normalized values.
#' It corresponds to the deviance residuals when the \code{W} is not included
#' in the model but the gene and cell-level covariates are. As a results, when
#' \code{W} is not included in the model, the deviance residuals should capture
#' the biology.
#'
#' @details If one has already fitted a model using \code{\link{ZinbModel}},
#' the object containing such model can be used as input of \code{zinbwave} to
#' save the resulting W into a \code{SummarizedExperiment} and optionally
#' compute residuals and normalized values, without the need for re-fitting the
#' model.
#'
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#'
#' @examples
#' se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
#'                            colData = data.frame(bio = gl(2, 3)))
#'
#' m <- zinbwave(se, X="~bio")
setMethod("zinbwave", "SummarizedExperiment",
          function(Y, X, V, fitted_model, commondispersion=TRUE, verbose=FALSE,
                   nb.repeat.initialize=2, maxiter.optimize=25,
                   stop.epsilon.optimize=.0001,
                   BPPARAM=BiocParallel::bpparam(),
                   normalizedValues = TRUE, residuals = FALSE,
                   imputedValues = FALSE, ...) {

              if(missing(fitted_model)) {
                  fitted_model <- zinbFit(Y, X, V, commondispersion,
                                          verbose, nb.repeat.initialize,
                                          maxiter.optimize,
                                          stop.epsilon.optimize, BPPARAM, ...)
              }

              out <- as(Y, "SingleCellExperiment")

              if (normalizedValues){
                  norm <- computeDevianceResiduals(fitted_model, t(assay(Y)),
                                                   ignoreW = TRUE)
                  assay(out, "normalizedValues") <- t(norm)
              }

              if (residuals){
                  devres <- computeDevianceResiduals(fitted_model, t(assay(Y)),
                                                     ignoreW = FALSE)
                  assay(out, "residuals") <- t(devres)
              }

              if (imputedValues){
                  imputed <- imputeZeros(fitted_model, t(assay(Y)))
                  assay(out, "imputedValues") <- t(imputed)
              }

              if (nFactors(fitted_model) > 0){
                  W <- getW(fitted_model)
                  colnames(W) <- paste0('W', seq_len(nFactors(fitted_model)))
                  reducedDim(out, "zinbwave") <- W
              }

              return(out)
          }
)
