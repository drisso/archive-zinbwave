#' Class ZinbModel
#' 
#' Objects of this class store all the values needed to work with a fitted 
#' model. Starting with an object of this class, one should be able to simulate 
#' data from a zero-inflated negative binomial and to compute the 
#' log-likelihood.
#' 
#' @slot X matrix. The design matrix containing sample-level covariates.
#' @slot V matrix. The design matrix containing gene-level covariates.
#' @slot O_mu matrix. The offset matrix for mu.
#' @slot O_pi matrix. The offset matrix for pi.
#' @slot which_X_mu integer. Indeces of which columns of X to use in the 
#'   regression of mu.
#' @slot which_V_mu integer. Indeces of which columns of V to use in the 
#'   regression of mu.
#' @slot which_X_pi integer. Indeces of which columns of X to use in the 
#'   regression of pi.
#' @slot which_V_pi integer. Indeces of which columns of V to use in the 
#'   regression of pi.
#' @slot W matrix. The factors of gene-level latent factors.
#' @slot beta_mu matrix or NULL. The coefficients of X in the regression of mu.
#' @slot gamma_mu matrix or NULL. The coefficients of V in the regression of mu.
#' @slot alpha_mu matrix or NULL. The coefficients of W in the regression of mu.
#' @slot beta_pi matrix or NULL. The coefficients of X in the regression of pi.
#' @slot gamma_pi matrix or NULL. The coefficients of V in the regression of pi.
#' @slot alpha_pi matrix or NULL. The coefficients of W in the regression of pi.
#' @slot logtheta numeric. A vector of log of inverse dispersion parameters.
#' @slot epsilon numeric. The regularization parameter for penalized maximum
#'   likelihood estimation
#' @slot penalty_alpha_mu numeric. This is a scalar that multiplies
#'   epsilon to allow differential regularization for the rows of alpha_mu
#'   parameters.
#' @slot penalty_beta_mu numeric. This is a vector that multiplies
#'   epsilon to allow differential regularization for the rows of beta_mu
#'   parameters.
#' @slot penalty_gamma_mu numeric. This is a vector that multiplies
#'   epsilon to allow differential regularization for the rows of gamma_mu
#'   parameters.
#' @slot penalty_alpha_pi numeric. This is a scalar that multiplies
#'   epsilon to allow differential regularization for the matrix of alpha_pi
#'   parameters.
#' @slot penalty_beta_pi numeric. This is a vector that multiplies
#'   epsilon to allow differential regularization for the rows of beta_pi
#'   parameters.
#' @slot penalty_gamma_pi numeric. This is a vector that multiplies
#'   epsilon to allow differential regularization for the rows of gamma_pi
#'   parameters.
#' @slot penalty_W numeric. This is a vector that multiplies epsilon to
#'   allow differential regularization for the columns of W parameters.
#' @slot epsilon_varphi numeric. This is the penalty for the variance of the
#'   dispersion parameter.
#'   
#' @details For the full description of the model see the model vignette. 
#'   Internally, the slots are checked so that the matrices are of the 
#'   appropriate dimensions: in particular, \code{X}, \code{O_mu}, \code{O_pi},
#'   and \code{W} need to have \code{n} rows, \code{V} needs to have \code{J}
#'   rows, \code{logtheta} must be of length \code{J}.
#' @name ZinbModel-class
#' @import methods
#' @exportClass ZinbModel
#' @docType Class
#' 
setClass(
    Class = "ZinbModel",
    slots = list(X = "matrix",
                 V = "matrix",
                 O_mu = "matrix",
                 O_pi = "matrix",
                 which_X_mu = "integer",
                 which_V_mu = "integer",
                 which_X_pi = "integer",
                 which_V_pi = "integer",
                 X_intercept = "logical",
                 V_intercept = "logical",
                 W = "matrix",
                 beta_mu = "matrix",
                 gamma_mu = "matrix",
                 alpha_mu = "matrix",
                 beta_pi = "matrix",
                 gamma_pi = "matrix",
                 alpha_pi = "matrix",
                 logtheta = "numeric",
                 epsilon = "numeric",
                 penalty_alpha_mu = "numeric",
                 penalty_beta_mu = "numeric",
                 penalty_gamma_mu = "numeric",
                 penalty_alpha_pi = "numeric",
                 penalty_beta_pi = "numeric",
                 penalty_gamma_pi = "numeric",
                 penalty_W = "numeric",
                 epsilon_varphi = "numeric"
                 )
)

setValidity("ZinbModel", function(object){
    n <- NROW(object@X) # number of samples
    J <- NROW(object@V) # number of genes
    K <- NCOL(object@W) # number of latent factors
    
    if(K > n) {
        return("Cannot have more latent factors than samples.")        
    }
    if(K > J) {
        return("Cannot have more latent factors than genes.")        
    }
    if(NROW(object@W) != n) {
        return("W must have n rows!")
    }
    if((length(object@which_X_mu)>0) && (max(object@which_X_mu) > NCOL(object@X))) {
        return("which_X_mu: subscript out of bound!")
    }
    if((length(object@which_X_pi)>0) && (max(object@which_X_pi) > NCOL(object@X))) {
        return("which_X_pi: subscript out of bound!")
    }
    if((length(object@which_V_mu)>0) && (max(object@which_V_mu) > NCOL(object@V))) {
        return("which_V_mu: subscript out of bound!")
    }
    if((length(object@which_V_pi)>0) && (max(object@which_V_pi) > NCOL(object@V))) {
        return("which_V_pi: subscript out of bound!")
    }
    if(length(object@which_X_mu) != NROW(object@beta_mu)){
        return("beta_mu must have the same number of rows as there are indices in which_X_mu!")
    }
    if(length(object@which_X_pi) != NROW(object@beta_pi)){
        return("beta_pi must have the same number of rows as there are indices in which_X_pi!")
    }
    if(length(object@which_V_mu) != NROW(object@gamma_mu)){
        return("gamma_mu must have the same number of rows as there are indices in which_V_mu!")
    }
    if(length(object@which_V_pi) != NROW(object@gamma_pi)){
        return("gamma_pi must have the same number of rows as there are indices in which_V_pi!")
    }
    if(NCOL(object@beta_mu) != J) {
        return("beta_mu must have J columns!")
    }
    if(NCOL(object@beta_pi) != J) {
        return("beta_pi must have J columns!")
    }
    if(NCOL(object@gamma_mu) != n) {
        return("gamma_mu must have n columns!")
    }
    if(NCOL(object@gamma_pi) != n) {
        return("gamma_pi must have n columns!")
    }
    if(NCOL(object@alpha_mu) != J) {
        return("alpha_mu must have J columns!")
    }
    if(NCOL(object@alpha_pi) != J) {
        return("alpha_pi must have J columns!")
    }
    if(NROW(object@alpha_mu) != K) {
        return("alpha_mu must have K rows!")
    }
    if(NROW(object@alpha_pi) != K) {
        return("alpha_pi must have K rows!")
    }
    if(NROW(object@O_mu) != n) {
        return("O_mu must have n rows!")
    }
    if(NROW(object@O_pi) != n) {
        return("O_pi must have n rows!")
    }
    if(NCOL(object@O_mu) != J) {
        return("O_mu must have J columns!")
    }
    if(NCOL(object@O_pi) != J) {
        return("O_pi must have J columns!")
    }
    if(length(object@logtheta) != J) {
        return("logtheta must have length J!")
    }
    if(length(object@penalty_alpha_mu) > 1) {
        return("penalty_alpha_mu must be a scalar !")
    }
    if(length(object@penalty_beta_mu) != NROW(object@beta_mu)) {
        return("The length of penalty_beta_mu must be the number of rows of beta_mu !")
    }
    if(length(object@penalty_gamma_mu) != NROW(object@gamma_mu)) {
        return("The length of penalty_gamma_mu must be the number of rows of gamma_mu !")
    }
    if(length(object@penalty_alpha_pi) > 1) {
        return("penalty_alpha_pi can only be a scalar !")
    }
    if(length(object@penalty_beta_pi) != NROW(object@beta_pi)) {
        return("The length of penalty_beta_pi must be the number of rows of beta_pi !")
    }
    if(length(object@penalty_gamma_pi) != NROW(object@gamma_pi)) {
        return("The length of penalty_gamma_pi must be the number of rows of gamma_pi !")
    }
    if(length(object@penalty_W) != NCOL(object@W)) {
        return("The length of penalty_W must be the number of columns of W !")
    }
    return(TRUE)
}
)