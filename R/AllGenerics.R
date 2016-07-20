#' Generic function that returns the number of latent factors
#' 
#' Given an object that describes a dataset or a model involving latent factors,
#' this function returns the number of latent factors.
#' @param x an object that describes a dataset or a model involving latent
#'   factors
#' @return the number of latent factors
#' @export
setGeneric("nFactors", function(x) standardGeneric("nFactors"))

#' Returns the matrix of mean parameters
#' 
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the matrix of mean parameters.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the matrix of mean parameters
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getMu(a) 
#' @export
setGeneric("getMu", function(object) standardGeneric("getMu"))

#' Returns the matrix of probabilities of zero
#' 
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the matrix of probabilities of 0.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the matrix of probabilities of 0
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getPi(a) 
#' @export
setGeneric("getPi", function(object) standardGeneric("getPi"))

#' Returns the vector of dispersion parameters
#' 
#' Given an object that describes a matrix of zero-inflated negative binomial
#' distributions, returns the vector of dispersion parameters \code{phi}.
#' @param object an object that describes a matrix of zero-inflated.
#'   distributions.
#' @return the vector of dispersion parameters
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getPhi(a) 
#' @export
setGeneric("getPhi", function(object) standardGeneric("getPhi"))

#' Returns the vector of inverse dispersion parameters
#' 
#' Given an object that describes a matrix of zero-inflated negative binomial 
#' distributions, returns the vector of inverse dispersion parameters
#' \code{theta}.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the vector of inverse dispersion parameters theta
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getTheta(a) 
#' @export
setGeneric("getTheta", function(object) standardGeneric("getTheta"))

#' Simulate counts from a zero-inflated negative binomial model
#' 
#' Given an object that describes zero-inflated negative binomial distribution,
#' simulate counts from the distribution.
#' @param object an object that describes a matrix of zero-inflated negative
#'   binomial.
#' @param seed an optional integer to specify how the random number generator
#'   should be initialized with a call to \code{set.seed}. If missing, the
#'   random generator state is not changed.
#' @param ... additional arguments.
#' @return A list with the following elements.
#'   \itemize{
#'   \item{counts}{the matrix with the simulated counts.}
#'   \item{dataNB}{the data simulated from the negative binomial.}
#'   \item{dataDropouts}{the data simulated from the binomial process.}
#'   \item{zeroFraction}{the fraction of zeros.}
#'   }
#' @examples
#' a <- zinbModel(n=5, J=10)
#' zinbSim(a) 
#' @export
setGeneric("zinbSim",function(object, seed, ...) standardGeneric("zinbSim"))

#' Compute the log-likelihood of a model given some data
#' 
#' Given a statistical model and some data, this function computes the
#' log-likelihood of the model given the data, i.e., the log-probability of the
#' data under the model.
#' @param model an object that describes a statistical model.
#' @param x an object that describes data.
#' @param ... additional arguments.
#' @return The log-likelihood of the model given the data.
#' @examples
#' m <- zinbModel(n=5, J=10)
#' x <- zinbSim(m)
#' loglik(m, x$counts)
#' @export
setGeneric("loglik", function(model, x, ...) standardGeneric("loglik"))

#' Compute the penalty of a model
#' 
#' Given a statistical model with regularization parameters, compute the
#' penalty.
#' @param model an object that describes a statistical model with regularization
#'   parameters.
#' @return The penalty of the model.
#' @examples
#' m <- zinbModel(K=2)
#' penalty(m)
#' @export
setGeneric("penalty", function(model) standardGeneric("penalty"))

#' Fit a ZINB regression model
#' 
#' Given an object with the data, it fits a ZINB model.
#' 
#' @param Y The data (genes in rows, samples in columns).
#' @param ... Additional parameters to describe the model.
#' @return An object of class \code{ZinbModel} that has been fitted by penalized
#'   maximum likelihood on the data.
#' @export
setGeneric("zinbFit", function(Y, ...) standardGeneric("zinbFit"))

#' Returns the sample-level design matrix for mu
#' 
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the sample-level design matrix for mu
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @param ... Additional parameters.
#' @return the sample-level design matrix for mu
#' @export
setGeneric("getX_mu", function(object, ...) standardGeneric("getX_mu"))

#' Returns the sample-level design matrix for pi
#' 
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the sample-level design matrix for pi
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @param ... Additional parameters.
#' @return the sample-level design matrix for pi
#' @export
setGeneric("getX_pi", function(object, ...) standardGeneric("getX_pi"))

#' Returns the gene-level design matrix for mu
#' 
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the gene-level design matrix for mu
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @param ... Additional parameters.
#' @return the gene-level design matrix for mu
#' @export
setGeneric("getV_mu", function(object, ...) standardGeneric("getV_mu"))

#' Returns the gene-level design matrix for pi
#' 
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the gene-level design matrix for pi
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @param ... Additional parameters.
#' @return the gene-level design matrix for pi
#' @export
setGeneric("getV_pi", function(object, ...) standardGeneric("getV_pi"))

#' Returns the regularization parameter for alpha_mu
#' 
#' The regularization parameter is given by epsilon multiplied by
#' penalty_alpha_mu
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @export
setGeneric("getEpsilon_alpha_mu", function(object) standardGeneric("getEpsilon_alpha_mu"))

#' Returns the regularization parameter for beta_mu
#' 
#' The regularization parameter is given by epsilon multiplied by
#' penalty_beta_mu
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @export
setGeneric("getEpsilon_beta_mu", function(object) standardGeneric("getEpsilon_beta_mu"))

#' Returns the regularization parameter for gamma_mu
#' 
#' The regularization parameter is given by epsilon multiplied by
#' penalty_gamma_mu
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @export
setGeneric("getEpsilon_gamma_mu", function(object) standardGeneric("getEpsilon_gamma_mu"))

#' Returns the regularization parameter for alpha_pi
#' 
#' The regularization parameter is given by epsilon multiplied by
#' penalty_alpha_pi
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @export
setGeneric("getEpsilon_alpha_pi", function(object) standardGeneric("getEpsilon_alpha_pi"))

#' Returns the regularization parameter for beta_pi
#' 
#' The regularization parameter is given by epsilon multiplied by
#' penalty_beta_pi
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @export
setGeneric("getEpsilon_beta_pi", function(object) standardGeneric("getEpsilon_beta_pi"))

#' Returns the regularization parameter for gamma_pi
#' 
#' The regularization parameter is given by epsilon multiplied by
#' penalty_gamma_pi
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @export
setGeneric("getEpsilon_gamma_pi", function(object) standardGeneric("getEpsilon_gamma_pi"))

#' Returns the regularization parameter for W
#' 
#' The regularization parameter is given by epsilon multiplied by
#' penalty_W
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @export
setGeneric("getEpsilon_W", function(object) standardGeneric("getEpsilon_W"))
