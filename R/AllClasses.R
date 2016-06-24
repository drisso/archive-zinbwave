setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))

#' Class zinb_model
#' 
#' Objects of this class store all the values needed to work with a fitted 
#' model. Starting with an object of this class, one should be able to simulate 
#' data from a zero-inflated negative binomial and to compute the
#' log-likelihood.
#' 
#' @slot n integer. The number of samples.
#' @slot J integer. The number of genes.
#' @slot k_mu integer. The number of factors of unwanted variation for mu.
#' @slot k_pi integer. The number of factors of unwanted variation for pi.
#' @slot X matrix. The design matrix containing sample-level covariates.
#' @slot V matrix. The design matrix containing gene-level covariates.
#' @slot O_mu matrix. The offsets for mu.
#' @slot O_pi matrix. The offsets for pi.
#' @slot which_X_mu integer. Indeces of which columns of X to use in the
#'   regression of mu.
#' @slot which_V_mu integer. Indeces of which columns of V to use in the
#'   regression of mu.
#' @slot which_X_pi integer. Indeces of which columns of X to use in the
#'   regression of pi.
#' @slot which_V_pi integer. Indeces of which columns of V to use in the
#'   regression of pi.
#' @slot W matrix. The factors of unwanted variation.
#' @slot beta_mu matrix or NULL. The coefficients of X in the regression of mu.
#' @slot gamma_mu matrix or NULL. The coefficients of V in the regression of mu.
#' @slot alpha_mu matrix or NULL. The coefficients of W in the regression of mu.
#' @slot beta_pi matrix or NULL. The coefficients of X in the regression of pi.
#' @slot gamma_pi matrix or NULL. The coefficients of V in the regression of pi.
#' @slot alpha_pi matrix or NULL. The coefficients of W in the regression of pi.
#' @slot phi numeric. A vector of dispersion parameters.
#' 
#' @details For the full description of the model see the model vignette.
#'   Internally, the slots are checked so that the matrices are of the
#'   appropriate dimensions: in particular, X, O_mu, O_pi, and W need to have n
#'   rows, V needs to have J rows, phi must be of length J.
#' @name zinb_model-class
#' @aliases zinb_model
#' @import methods
#' @export
#' 
setClass(
    Class = "zinb_model",
    slots = list(n = "integer",
                 J = "integer",
                 k_mu = "integer",
                 k_pi = "integer",
                 X = "matrix",
                 V = "matrix",
                 O_mu = "matrix",
                 O_pi = "matrix",
                 which_X_mu = "integer",
                 which_V_mu = "integer",
                 which_X_pi = "integer",
                 which_V_pi = "integer",
                 W = "matrix",
                 beta_mu = "matrixOrNULL",
                 gamma_mu = "matrixOrNULL",
                 alpha_mu = "matrixOrNULL",
                 beta_pi = "matrixOrNULL",
                 gamma_pi = "matrixOrNULL",
                 alpha_pi = "matrixOrNULL",
                 phi = "numeric"
                 )
)

setValidity("zinb_model", function(object){
    if(k_mu > n | k_pi > n) {
        return("Cannot have more factors of unwanted variation than samples.")        
    }
    
    if(NROW(X) != n) {
        return("X must have n rows!")
    }
    if(NROW(V) != J) {
        return("V must have J rows!")
    }
    if(NROW(W) > 0 & NROW(W) != n) {
        return("If specified, W must have n rows!")
    }
    
    if(max(which_X_mu) > NCOL(X)) {
        return("which_X_mu: subscript out of bound!")
    }
    if(max(which_X_pi) > NCOL(X)) {
        return("which_X_pi: subscript out of bound!")
    }
    if(max(which_V_mu) > NCOL(V)) {
        return("which_V_mu: subscript out of bound!")
    }
    if(max(which_V_pi) > NCOL(V)) {
        return("which_V_pi: subscript out of bound!")
    }
    
    return(TRUE)
}
)