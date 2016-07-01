setMethod(
    f="initialize",
    signature="ZinbModel",
    definition=function(.Object, X, V, O_mu, O_pi, which_X_mu,
                        which_X_pi, which_V_mu, which_V_pi, W,beta_mu, beta_pi,
                        gamma_mu, gamma_pi, alpha_mu, alpha_pi, phi,epsilon,
                        penalty_alpha_mu, penalty_beta_mu, penalty_gamma_mu,
                        penalty_alpha_pi, penalty_beta_pi, penalty_gamma_pi,
                        penalty_W, epsilon_varphi, n, J, K) {
        
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
            } else if (!missing(phi)) {
                n <- length(phi)
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
        
        # Set the different slots
        if (!missing(X)) {
            .Object@X <- X
        } else {
            .Object@X <- matrix(nrow=n, ncol=0)
        }
        if (!missing(V)) {
            .Object@V <- V
        } else {
            .Object@V <- matrix(nrow=J, ncol=0)
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
        if (!missing(phi)) {
            .Object@logtheta <- -log(phi)
        } else {
            .Object@logtheta <- numeric(J)
        }
        if (!missing(epsilon)) {
            .Object@epsilon <- epsilon
        } else {
            .Object@epsilon <- 1e-3
        }
        if (!missing(penalty_alpha_mu)) {
            .Object@penalty_alpha_mu <- penalty_alpha_mu
        } else {
            val <- max(1, NROW(.Object@alpha_mu)*NCOL(.Object@alpha_mu))
            .Object@penalty_alpha_mu <- rep(1/val,NROW(.Object@alpha_mu))
        }
        if (!missing(penalty_beta_mu)) {
            .Object@penalty_beta_mu <- penalty_beta_mu
        } else {
            val <- max(1, NROW(.Object@beta_mu)*NCOL(.Object@beta_mu))
            .Object@penalty_beta_mu <- rep(1/val,NROW(.Object@beta_mu))
        }
        if (!missing(penalty_gamma_mu)) {
            .Object@penalty_gamma_mu <- penalty_gamma_mu
        } else {
            val <- max(1, NROW(.Object@gamma_mu)*NCOL(.Object@gamma_mu))
            .Object@penalty_gamma_mu <- rep(1/val,NROW(.Object@gamma_mu))
        }
        if (!missing(penalty_alpha_pi)) {
            .Object@penalty_alpha_pi <- penalty_alpha_pi
        } else {
            val <- max(1, NROW(.Object@alpha_pi)*NCOL(.Object@alpha_pi))
            .Object@penalty_alpha_pi <- rep(1/val,NROW(.Object@alpha_pi))
        }
        if (!missing(penalty_beta_pi)) {
            .Object@penalty_beta_pi <- penalty_beta_pi
        } else {
            val <- max(1, NROW(.Object@beta_pi)*NCOL(.Object@beta_pi))
            .Object@penalty_beta_pi <- rep(1/val,NROW(.Object@beta_pi))
        }
        if (!missing(penalty_gamma_pi)) {
            .Object@penalty_gamma_pi <- penalty_gamma_pi
        } else {
            val <- max(1, NROW(.Object@gamma_pi)*NCOL(.Object@gamma_pi))
            .Object@penalty_gamma_pi <- rep(1/val,NROW(.Object@gamma_pi))
        }
        if (!missing(penalty_W)) {
            .Object@penalty_W <- penalty_W
        } else {
            val <- max(1, NROW(.Object@W)*NCOL(.Object@W))
            .Object@penalty_W <- rep(1/val,NCOL(.Object@W))
        }
        if (!missing(epsilon_varphi)) {
            .Object@epsilon_varphi <- epsilon_varphi
        } else {
            .Object@epsilon_varphi <- 1
        }
        
        validObject(.Object) # call of the inspector
        return(.Object)
    }
)

#' Initialize an object of class ZinbModel
#' @export
#' 
#' @param ... arguments passed to \code{new()}. See the \code{slots} section in 
#'   \code{\link{`ZinbModel-class`}}.
#'   
#' @details This is a light wrapper around the new() function to create an 
#'   instance of class \code{ZinbModel}.
#'   
#' @details If any of the related matrices are passed, \code{n}, \code{J}, and
#'   \code{K} are inferred. Alternatively, the user can specify one or more of
#'   \code{n}, \code{J}, and \code{K}.
#'   
#' @details A call with no argument has the following default values: \code{n =
#'   50}, \code{J = 100}, \code{K = 2}.
#'   
#' @details Although it is possible to create new instances of the class by 
#'   calling this function, this is not the most common way of creating
#'   \code{ZinbModel} objects. The main use of the class is within the
#'   \code{\link{zinbFit}} function.
#' 
#' @examples
#' a <- zinbModel()
#' getN(a)
#' getJ(a)
#' getK(a)
#' 
zinbModel <- function(...) {
    new(Class="ZinbModel", ...)
}

#' Returns the number of samples
#' 
#' Given an object that describes a ZINB model, this function returns the number
#' of samples.
#' @param x an object of class \code{ZinbModel}.
#' @return the number of samples
#' @export
#' @importFrom clusterExperiment nSamples
#' @examples
#' a <- zinbModel()
#' nSamples(a) 
setMethod("nSamples", "ZinbModel",
          function(x) {
              return(NROW(x@X))
          }
)

#' Returns the number of features
#' 
#' Given an object that describes a ZINB model, this
#' function returns the number of genes.
#' @param x an object of class \code{ZinbModel}.
#' @return the number of genes.
#' @examples
#' a <- zinbModel()
#' nFeatures(a)
#' @export
#' @importFrom clusterExperiment nFeatures
setMethod("nFeatures", "ZinbModel",
          function(x) {
              return(NROW(x@V))
          }
)

#' @export
#' @describeIn getK return the number of latent factors.
setMethod("getK", "ZinbModel",
          function(object) {
              return(NCOL(object@W))
          }
)

#' @export
#' @describeIn getMu return the mean of the non-zero component.
setMethod("getMu", "ZinbModel",
    function(object) {
        return(exp(object@X[,object@which_X_mu] %*% object@beta_mu +
                       t(object@V[,object@which_V_mu] %*% object@gamma_mu) + 
                       object@O_mu))
    }
)

#' @export
#' @describeIn getPi return the probability of zero.
#' @importFrom stats binomial
setMethod("getPi", "ZinbModel",
    function(object) {
      return(stats::binomial()$linkinv(object@X[,object@which_X_mu] %*% object@beta_mu + 
               t(object@V[,object@which_V_mu] %*% object@gamma_mu) + object@O_mu))
    }
)

#' @export
#' @describeIn getPhi return the dispersion parameter.
setMethod("getPhi", "ZinbModel",
          function(object) {
              return(exp(-object@logtheta))
          }
)

#' @export
#' @describeIn getTheta return the inverse of the dispersion parameter.
setMethod("getTheta", "ZinbModel",
          function(object) {
              return(exp(object@logtheta))
          }
)

#' @export
#' @describeIn simulateZINB simulate from a ZINB distribution.
#' @importFrom parallel mclapply
setMethod(
    f="simulateZINB",
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
        n <- getN(object)
        J <- getJ(object)
        
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
        zinb.loglik(x, getMu(model), getTheta(model), getPi(model))
    }
)

#' @export
#' @describeIn penalty return the penalization.
setMethod(
    f="penalty",
    signature="ZinbModel",
    definition=function(model) {
        sum(model@epsilon*model@penalty_alpha_mu*(model@alpha_mu)^2)/2
        + sum(model@epsilon*model@penalty_alpha_pi*(model@alpha_pi)^2)/2
        + sum(model@epsilon*model@penalty_beta_mu*(model@beta_mu)^2)/2
        + sum(model@epsilon*model@penalty_beta_pi*(model@beta_pi)^2)/2
        + sum(model@epsilon*model@penalty_gamma_mu*(model@gamma_mu)^2)/2
        + sum(model@epsilon*model@penalty_gamma_pi*(model@gamma_pi)^2)/2
        + sum(model@epsilon*model@penalty_W*(model@W)^2)/2
        + model@epsilon_varphi*var(getPhi(model))
    }
)