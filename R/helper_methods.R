#' Toy dataset to check the model
#' 
#' @name toydata
#' @aliases toydata
#' 
#' @format A matrix of integers (counts) with 96 samples (rows) and 500 genes
#'   (columns).
NULL

setMethod(
    f="initialize",
    signature="ZinbModel",
    definition=function(.Object, X, V, O_mu, O_pi, which_X_mu,
                        which_X_pi, which_V_mu, which_V_pi, W, beta_mu, beta_pi,
                        gamma_mu, gamma_pi, alpha_mu, alpha_pi, zeta, epsilon,
                        epsilon_beta_mu, epsilon_gamma_mu, epsilon_beta_pi, 
                        epsilon_gamma_pi, epsilon_W, epsilon_alpha, 
                        epsilon_zeta, epsilon_min_logit, n, J, K) {
        
        # Find n (default 50), J (default 100), K (default 0)
        if (missing(n)) {
            if (!missing(X)) {
                n <- NROW(X)
            } else if (!missing(gamma_mu)) {
                n <- NCOL(gamma_mu)
            } else if (!missing(gamma_pi)) {
                n <- NCOL(gamma_pi)
            } else if (!missing(W)) {
                n <- NROW(W)
            } else if (!missing(O_mu)) {
                n <- NROW(O_mu)
            } else if (!missing(O_pi)) {
                n <- NROW(O_pi)
            } else {
                n <- 50
            }
        }
        if (missing(J)) {
            if (!missing(V)) {
                J <- NROW(V)
            } else if (!missing(beta_mu)) {
                J <- NCOL(beta_mu)
            } else if (!missing(beta_pi)) {
                J <- NCOL(beta_pi)
            } else if (!missing(alpha_mu)) {
                J <- NCOL(alpha_mu)
            } else if (!missing(alpha_pi)) {
                J <- NCOL(alpha_pi)
            } else if (!missing(O_mu)) {
                J <- NCOL(O_mu)
            } else if (!missing(O_pi)) {
                J <- NCOL(O_pi)
            } else if (!missing(zeta)) {
                J <- length(zeta)
            } else {
                J <- 100
            }
        }
        if (missing(K)) {
            if (!missing(W)) {
                K <- NCOL(W)
            } else {
                K <- 0
            }
        }
        
        # Set the different slots for the matrices
        if (!missing(X)) {
            .Object@X <- X
        } else {
            .Object@X <- matrix(1, nrow=n, ncol=1)
        }
        if (!missing(V)) {
            .Object@V <- V
        } else {
            .Object@V <- matrix(1, nrow=J, ncol=1)
        }
        
        if (!missing(which_X_mu)) {
            .Object@which_X_mu <- which_X_mu
        } else if (NCOL(.Object@X)>0) {
            .Object@which_X_mu <- seq(NCOL(.Object@X))
        } else {
            .Object@which_X_mu <- integer(0)
        }
        if (!missing(which_X_pi)) {
            .Object@which_X_pi <- which_X_pi
        } else if (NCOL(.Object@X)>0) {
            .Object@which_X_pi <- seq(NCOL(.Object@X))
        } else {
            .Object@which_X_pi <- integer(0)
        }
        if (!missing(which_V_mu)) {
            .Object@which_V_mu <- which_V_mu
        } else if (NCOL(.Object@V)>0) {
            .Object@which_V_mu <- seq(NCOL(.Object@V))
        } else {
            .Object@which_V_mu <- integer(0)
        }
        if (!missing(which_V_pi)) {
            .Object@which_V_pi <- which_V_pi
        } else if (NCOL(.Object@V)>0) {
            .Object@which_V_pi <- seq(NCOL(.Object@V))
        } else {
            .Object@which_V_pi <- integer(0)
        }
        
        
        if (all(.Object@X[,.Object@which_X_mu[1]]==1)) {
            .Object@X_mu_intercept <- TRUE
        } else {
            .Object@X_mu_intercept <- FALSE
        }
        
        if (all(.Object@V[,.Object@which_V_mu[1]]==1)) {
            .Object@V_mu_intercept <- TRUE
        } else {
            .Object@V_mu_intercept <- FALSE
        }
        if (all(.Object@X[,.Object@which_X_pi[1]]==1)) {
            .Object@X_pi_intercept <- TRUE
        } else {
            .Object@X_pi_intercept <- FALSE
        }
        
        if (all(.Object@V[,.Object@which_V_pi[1]]==1)) {
            .Object@V_pi_intercept <- TRUE
        } else {
            .Object@V_pi_intercept <- FALSE
        }
        
        if (!missing(O_mu)) {
            .Object@O_mu <- O_mu
        } else {
            .Object@O_mu <- matrix(0, nrow=n, ncol=J)
        }
        if (!missing(O_pi)) {
            .Object@O_pi <- O_pi
        } else {
            .Object@O_pi <- matrix(0, nrow=n, ncol=J)
        }
        
        if (!missing(W)) {
            .Object@W <- W
        } else {
            .Object@W <- matrix(0, nrow=n , ncol=K)
        }
        
        if (!missing(beta_mu)) {
            .Object@beta_mu <- beta_mu
        } else {
            .Object@beta_mu <- matrix(0, nrow=length(.Object@which_X_mu) , ncol=J)
        }
        if (!missing(beta_pi)) {
            .Object@beta_pi <- beta_pi
        } else {
            .Object@beta_pi <- matrix(0, nrow=length(.Object@which_X_pi) , ncol=J)
        }
        if (!missing(gamma_mu)) {
            .Object@gamma_mu <- gamma_mu
        } else {
            .Object@gamma_mu <- matrix(0, nrow=length(.Object@which_V_mu) , ncol=n)
        }
        if (!missing(gamma_pi)) {
            .Object@gamma_pi <- gamma_pi
        } else {
            .Object@gamma_pi <- matrix(0, nrow=length(.Object@which_V_pi) , ncol=n)
        }
        if (!missing(alpha_mu)) {
            .Object@alpha_mu <- alpha_mu
        } else {
            .Object@alpha_mu <- matrix(0, nrow=K , ncol=J)
        }
        if (!missing(alpha_pi)) {
            .Object@alpha_pi <- alpha_pi
        } else {
            .Object@alpha_pi <- matrix(0, nrow=K , ncol=J)
        }
        if (!missing(zeta)) {
            .Object@zeta <- zeta
        } else {
            .Object@zeta <- numeric(J)
        }
        
        # Regularization parameters
        if (missing(epsilon)) {
            epsilon <- 1e-3
        }
        if (missing(epsilon_min_logit)) {
            .Object@epsilon_min_logit <- 1e-3
        }
        if (!missing(epsilon_beta_mu)) {
            .Object@epsilon_beta_mu <- epsilon_beta_mu
        } else {
            .Object@epsilon_beta_mu <- epsilon/J
        }
        if (!missing(epsilon_gamma_mu)) {
            .Object@epsilon_gamma_mu <- epsilon_gamma_mu
        } else {
            .Object@epsilon_gamma_mu <- epsilon/n
        }
        if (!missing(epsilon_beta_pi)) {
            .Object@epsilon_beta_pi <- epsilon_beta_pi
        } else {
            .Object@epsilon_beta_pi <- epsilon/J
        }
        if (!missing(epsilon_gamma_pi)) {
            .Object@epsilon_gamma_pi <- epsilon_gamma_pi
        } else {
            .Object@epsilon_gamma_pi <- epsilon/n
        }        
        if (!missing(epsilon_W)) {
            .Object@epsilon_W <- epsilon_W
        } else {
            .Object@epsilon_W <- epsilon/n
        }        
        if (!missing(epsilon_alpha)) {
            .Object@epsilon_alpha <- epsilon_alpha
        } else {
            .Object@epsilon_alpha <- epsilon/J
        }
        if (!missing(epsilon_zeta)) {
            .Object@epsilon_zeta <- epsilon_zeta
        } else {
            .Object@epsilon_zeta <- epsilon
        }        
        
        validObject(.Object) # call of the inspector
        return(.Object)
    }
)

#' Initialize an object of class ZinbModel
#' @export
#' 
#' @param ... arguments passed to \code{new()}. See the \code{slots} section in 
#'   \code{\link{ZinbModel}}.
#'   
#' @details This is a light wrapper around the new() function to create an 
#'   instance of class \code{ZinbModel}.
#'   
#' @details If any of the related matrices are passed, \code{n}, \code{J}, and 
#'   \code{K} are inferred. Alternatively, the user can specify one or more of 
#'   \code{n}, \code{J}, and \code{K}.
#'   
#' @details The regularization parameters can be set by a unique parameter
#'   \code{epsilon}, as described in the vignette; or specific values for the
#'   different regularization parameters can also be provided.
#'   
#' @details A call with no argument has the following default values: \code{n = 
#'   50}, \code{J = 100}, \code{K = 0}, \code{epsilon=1e-3}.
#'   
#' @details Although it is possible to create new instances of the class by 
#'   calling this function, this is not the most common way of creating 
#'   \code{ZinbModel} objects. The main use of the class is within the 
#'   \code{\link{zinbFit}} function.
#'   
#' @examples
#' a <- zinbModel()
#' nSamples(a)
#' nFeatures(a)
#' nFactors(a)
#' 
zinbModel <- function(...) {
    new(Class="ZinbModel", ...)
}

#' @export
#' @describeIn ZinbModel show useful info on the object.
#' 
#' @param object an object of class \code{ZinbModel}.
setMethod("show", "ZinbModel",
          function(object) {
              cat(paste0("Object of class ZinbModel.\n",
                         NROW(object@X), " samples; ", NROW(object@V), " genes.\n",
                         NCOL(object@X), " sample-level covariates; ",
                         NCOL(object@V), " gene-level covariates; ",
                         NCOL(object@W), " latent factors.\n"))
          }
)


################################################################
# Extract various informations and variables from a ZINB model #
################################################################

#' @export
#' @importFrom clusterExperiment nSamples
#' @describeIn ZinbModel returns the number of samples.
#' @param x an object of class \code{ZinbModel}.
setMethod("nSamples", "ZinbModel",
          function(x) {
              return(NROW(x@X))
          }
)

#' @export
#' @importFrom clusterExperiment nFeatures
#' @describeIn ZinbModel returns the number of features.
setMethod("nFeatures", "ZinbModel",
          function(x) {
              return(NROW(x@V))
          }
)

#' @export
#' @describeIn ZinbModel returns the number of latent factors.
setMethod("nFactors", "ZinbModel",
          function(x) {
              return(NCOL(x@W))
          }
)


#' @export
#' @describeIn getX_mu return the sample-level design matrix for mu.
#' @param intercept logical. Whether to return the intercept (ignored if X_mu has
#'   no intercept). Default \code{TRUE}
setMethod("getX_mu", "ZinbModel",
          function(object, intercept=TRUE) {
              if(object@X_mu_intercept && !intercept) {
                  which_X_mu <- object@which_X_mu[-1]
              } else {
                  which_X_mu <- object@which_X_mu
              }
              return(object@X[, which_X_mu, drop=FALSE])
          }
)

#' @export
#' @describeIn getX_pi return the sample-level design matrix for pi.
#' @param intercept logical. Whether to return the intercept (ignored if X_pi has
#'   no intercept).
setMethod("getX_pi", "ZinbModel",
          function(object, intercept=TRUE) {
              if(object@X_pi_intercept && !intercept) {
                  which_X_pi <- object@which_X_pi[-1]
              } else {
                  which_X_pi <- object@which_X_pi
              }
              return(object@X[, which_X_pi, drop=FALSE])
          }
)

#' @export
#' @describeIn getV_mu return the gene-level design matrix for mu.
#' @param intercept logical. Whether to return the intercept (ignored if V_mu has
#'   no intercept).
setMethod("getV_mu", "ZinbModel",
          function(object, intercept=TRUE) {
              if(object@V_mu_intercept && !intercept) {
                  which_V_mu <- object@which_V_mu[-1]
              } else {
                  which_V_mu <- object@which_V_mu
              }
              return(object@V[, which_V_mu, drop=FALSE])
          }
)

#' @export
#' @describeIn getV_pi return the sample-level design matrix for pi.
#' @param intercept logical. Whether to return the intercept (ignored if V_pi has
#'   no intercept).
setMethod("getV_pi", "ZinbModel",
          function(object, intercept=TRUE) {
              if(object@V_pi_intercept && !intercept) {
                  which_V_pi <- object@which_V_pi[-1]
              } else {
                  which_V_pi <- object@which_V_pi
              }
              return(object@V[, which_V_pi, drop=FALSE])
          }
)

#' @export
#' @describeIn getLogMu return the logarithm of the mean of the non-zero
#'   component.
setMethod("getLogMu", "ZinbModel",
          function(object) {
              return(getX_mu(object) %*% object@beta_mu +
                         t(getV_mu(object) %*% object@gamma_mu) +
                         object@W %*% object@alpha_mu +
                         object@O_mu)
          }
)

#' @export
#' @describeIn getMu return the mean of the non-zero component.
setMethod("getMu", "ZinbModel",
    function(object) {
        return(exp(getLogMu(object)))
    }
)

#' @export
#' @describeIn getLogitPi return the logit-probability of zero.
#' @importFrom stats binomial
setMethod("getLogitPi", "ZinbModel",
          function(object) {
              return(getX_pi(object) %*% object@beta_pi +
                         t(getV_pi(object) %*% object@gamma_pi) +
                         object@W %*% object@alpha_pi +
                         object@O_pi)
          }
)

#' @export
#' @describeIn getPi return the probability of zero.
#' @importFrom stats binomial
setMethod("getPi", "ZinbModel",
    function(object) {
        # return(stats::binomial()$linkinv(getLogitPi(object))
        # Instead of the call to stats::binomial() in the previous line, we 
        # directly compute with the exp() function which remains exact for 
        # smaller values of the arguments.
        return(1/(1+exp(-getLogitPi(object))))
    }
)

#' @export
#' @describeIn getZeta return the log of the inverse of the dispersion parameter.
setMethod("getZeta", "ZinbModel",
          function(object) {
              return(object@zeta)
          }
)

#' @export
#' @describeIn getPhi return the dispersion parameter.
setMethod("getPhi", "ZinbModel",
          function(object) {
              return(exp(-object@zeta))
          }
)

#' @export
#' @describeIn getTheta return the inverse of the dispersion parameter.
setMethod("getTheta", "ZinbModel",
          function(object) {
              return(exp(object@zeta))
          }
)

#' @export
#' @describeIn getEpsilon_beta_mu method for ZinbModel.
setMethod("getEpsilon_beta_mu", "ZinbModel",
          function(object) {
              e <- rep(object@epsilon_beta_mu, length(object@which_X_mu))
              if (object@X_mu_intercept) {
                  e[1] <- 0
              }
              e
          }
)

#' @export
#' @describeIn getEpsilon_gamma_mu method for ZinbModel.
setMethod("getEpsilon_gamma_mu", "ZinbModel",
          function(object) {
              e <- rep(object@epsilon_gamma_mu, length(object@which_V_mu))
              if (object@V_mu_intercept) {
                  e[1] <- 0
              }
              e
          }
)

#' @export
#' @describeIn getEpsilon_beta_pi method for ZinbModel.
setMethod("getEpsilon_beta_pi", "ZinbModel",
          function(object) {
              e <- rep(object@epsilon_beta_pi, length(object@which_X_pi))
              if (object@X_pi_intercept) {
                  e[1] <- object@epsilon_min_logit
              }
              e
          }
)

#' @export
#' @describeIn getEpsilon_gamma_pi method for ZinbModel.
setMethod("getEpsilon_gamma_pi", "ZinbModel",
          function(object) {
              e <- rep(object@epsilon_gamma_pi, length(object@which_V_pi))
              if (object@V_pi_intercept) {
                  e[1] <- object@epsilon_min_logit
              }
              e
          }
)

#' @export
#' @describeIn getEpsilon_W method for ZinbModel.
setMethod("getEpsilon_W", "ZinbModel",
          function(object) {
              rep(object@epsilon_W, nFactors(object))
          }
)

#' @export
#' @describeIn getEpsilon_alpha method for ZinbModel.
setMethod("getEpsilon_alpha", "ZinbModel",
          function(object) {
              rep(object@epsilon_alpha, nFactors(object))
          }
)

#' @export
#' @describeIn getEpsilon_zeta method for ZinbModel.
setMethod("getEpsilon_zeta", "ZinbModel",
          function(object) {
              object@epsilon_zeta
          }
)

########################
# Other useful methods #
########################

#' @export
#' @describeIn zinbSim simulate from a ZINB distribution.
#' @importFrom parallel mclapply
#' 
#' @param no_cores number of cores for parallel computations (to be passed to
#'   mclapply).
setMethod(
    f="zinbSim",
    signature="ZinbModel",
    definition=function(object, seed, no_cores=1) {
        
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            runif(1)
        }
        
        if (missing(seed)) {
            RNGstate <- get(".Random.seed", envir = .GlobalEnv)
        } else {
            R.seed <- get(".Random.seed", envir = .GlobalEnv)
            set.seed(seed)
            RNGstate <- structure(seed, kind = as.list(RNGkind()))
            on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        }
        
        mu <- getMu(object)
        pi <- getPi(object)
        theta <- getTheta(object)
        n <- nSamples(object)
        J <- nFeatures(object)
        
        # Simulate negative binomial with the mean matrix and dispersion
        # parameters
        datanb <- parallel::mclapply(seq(n*J), 
                                     function(i) { 
                                         rnbinom(1, mu = mu[i] , size = theta[ceiling(i/n)])
                                     }, mc.cores = no_cores)
        
        data.nb <- matrix(unlist(datanb), nrow = n )
        
        # Simulate the binary dropout matrix. "1" means that a dropout (zero) is
        # observed instead of the value
        datado <- parallel::mclapply(seq(n*J), 
                                     function(i) { 
                                         rbinom(1, size =1, prob = pi[i])
                                     }, mc.cores = no_cores)
        
        data.dropout <- matrix(unlist(datado), nrow = n)
        
        # Matrix of zero-inflated counts
        counts <- data.nb * (1 - data.dropout)
        
        # Fraction of zeros in the matrix
        zero.fraction <- sum(counts == 0) / (n*J)    
        
        ret <- list(counts = counts, dataNB = data.nb, 
                    dataDropouts = data.dropout, zeroFraction = zero.fraction)
        attr(ret, "seed") <- RNGstate
        ret
    }
)

#' @export
#' @describeIn loglik return the log-likelihood of the ZINB model.
setMethod(
    f="loglik",
    signature=c("ZinbModel","matrix"),
    definition=function(model, x) {
        zinb.loglik(x, getMu(model), rep(getTheta(model),rep(nSamples(model),nFeatures(model))), getLogitPi(model))
    }
)

#' @export
#' @describeIn penalty return the penalization.
setMethod(
    f="penalty",
    signature="ZinbModel",
    definition=function(model) {
        sum(getEpsilon_alpha(model)*(model@alpha_mu)^2)/2 +
        sum(getEpsilon_alpha(model)*(model@alpha_pi)^2)/2 +
        sum(getEpsilon_beta_mu(model)*(model@beta_mu)^2)/2 +
        sum(getEpsilon_beta_pi(model)*(model@beta_pi)^2)/2 +
        sum(getEpsilon_gamma_mu(model)*(model@gamma_mu)^2)/2 +
        sum(getEpsilon_gamma_pi(model)*(model@gamma_pi)^2)/2 +
        sum(getEpsilon_W(model)*t(model@W)^2)/2 +
        getEpsilon_zeta(model)*var(getZeta(model))/2
    }
)

