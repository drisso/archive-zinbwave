#' Returns the number of samples
#' 
#' Given an object that describes a dataset or a model involving samples, this
#' function returns the number of samples.
#' @param object an object that describes a dataset or a model involving samples.
#' @return the number of samples
#' @examples
#' a <- zinb_model()
#' getN(a) 
#' @export
setGeneric("getN", function(object) standardGeneric("getN"))

#' Returns the number of genes
#' 
#' Given an object that describes a dataset or a model involving genes, this
#' function returns the number of genes.
#' @param object an object that describes a dataset or a model involving genes.
#' @return the number of genes.
#' @examples
#' a <- zinb_model()
#' getJ(a) 
#' @export
setGeneric("getJ", function(object) standardGeneric("getJ"))

#' Returns the number of latent factors
#' 
#' Given an object that describes a dataset or a model involving latent factors,
#' this function returns the number of latent factors.
#' @param object an object that describes a dataset or a model involving latent
#'   factors
#' @return the number of latent factors
#' @examples
#' a <- zinb_model()
#' getK(a) 
#' @export
setGeneric("getK", function(object) standardGeneric("getK"))

#' Returns the matrix of mean parameters
#' 
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the matrix of mean parameters.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the matrix of mean parameters
#' @examples
#' a <- zinb_model(n=5, J=10)
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
#' a <- zinb_model(n=5, J=10)
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
#' a <- zinb_model(n=5, J=10)
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
#' a <- zinb_model(n=5, J=10)
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
#' @return the matrix of probabilities of 0.
#' @examples
#' a <- zinb_model(n=5, J=10)
#' simulateZINB(a) 
#' @export
setGeneric("simulateZINB",function(object, seed, ...) standardGeneric("simulateZINB"))

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
#' m <- zinb_model(n=5, J=10)
#' x <- simulateZINB(m)
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
#' m <- zinb_model(K=2)
#' penalty(m)
#' @export
setGeneric("penalty", function(model) standardGeneric("penalty"))