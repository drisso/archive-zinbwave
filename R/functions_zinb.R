#' Fit a ZINB regression model
#'
#' @param Y The matrix of counts.
#' @param ... Additional parameters to describe the model
#' @return An object of class zinb_model that has been fitted by penalized maximum likelihood on the count matrix.
#' @examples 
#' m=zinbFit(matrix(10,3,5))
#' @export
zinbFit <- function(Y, ...) {
    
    n = NROW(Y) # number of samples
    J = NCOL(Y) # number of genes
    
    # Create a zinb_model object
    m = zinb_model(n=n, J=J, ...)
    
    # Initialize the parameters
    m = zinb_initialize(m,Y)
    
    # Optimize parameters
    m = zinb_optimize(m,Y)
    
    m
}

#' Initialize the parameters of a ZINB regression model
#'
#' The initialization performs quick optimization of the parameters with several simplifying assumptions compared to the true model: non-zero counts are models as log-Gaussian, zeros are modeled as dropouts. The dispersion parameter is not modified.
#' @param m The model of class zinb_model
#' @param Y The matrix of counts.
#' @return An object of class zinb_model similar to the one given as argument with modified parameters alpha_mu, alpha_pi, beta_mu, beta_pi, gamma_mu, gamma_pi, W.
#' @examples
#' Y = matrix(10,3,5)
#' m = zinb_model(n=NROW(Y),J=NCOL(Y))
#' m = zinb_initialize(m,Y)
#' @export
zinb_initialize <- function(m, Y) {

    # TODO: write this function    
    m
}

#' Optimize the parameters of a ZINB regression model
#'
#' The parameters of the model given as argument are optimized by penalized maximum likelihood on the count matrix given as argument. It is recommended to call zinb_initialize before this function to have good starting point for optimization, since the optimization problem is not convex and can only converge to a local minimum.
#' @param m The model of class zinb_model
#' @param Y The matrix of counts.
#' @return An object of class zinb_model similar to the one given as argument with modified parameters alpha_mu, alpha_pi, beta_mu, beta_pi, gamma_mu, gamma_pi, W.
#' @examples
#' Y = matrix(10,3,5)
#' m = zinb_model(n=NROW(Y),J=NCOL(Y))
#' m = zinb_initialize(m,Y)
#' m = zinb_optimize(m,Y)
#' @export
zinb_optimize <- function(m, Y) {
    
    # TODO: write this function    
    m
}
