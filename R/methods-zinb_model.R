setMethod(
    f="initialize",
    signature="zinb_model",
    definition=function(.Object,X,V,O_mu,O_pi,which_X_mu,which_X_pi,which_V_mu,which_V_pi,W,beta_mu,beta_pi,gamma_mu,gamma_pi,alpha_mu,alpha_pi,phi,n,J,K) {
        # Find n (default 50), J (default 100), K (default 0)
        if (missing(n)) {
            if (!missing(X)) {
                n = NROW(X)
            } else if (!missing(gamma_mu)) {
                n = NCOL(gamma_mu)
            } else if (!missing(gamma_pi)) {
                n = NCOL(gamma_pi)
            } else if (!missing(W)) {
                n = NROW(W)
            } else if (!missing(O_mu)) {
                n = NROW(O_mu)
            } else if (!missing(O_pi)) {
                n = NROW(O_pi)
            } else if (!missing(phi)) {
                n = length(phi)
            } else {
                n = 50
            }
        }
        if (missing(J)) {
            if (!missing(V)) {
                J = NROW(V)
            } else if (!missing(beta_mu)) {
                J = NCOL(beta_mu)
            } else if (!missing(beta_pi)) {
                J = NCOL(beta_pi)
            } else if (!missing(alpha_mu)) {
                J = NCOL(alpha_mu)
            } else if (!missing(alpha_pi)) {
                J = NCOL(alpha_pi)
            } else if (!missing(O_mu)) {
                J = NCOL(O_mu)
            } else if (!missing(O_pi)) {
                J = NCOL(O_pi)
            } else {
                J = 100
            }
        }    
        if (missing(K)) {
            if (!missing(W)) {
                K = NCOL(W)
            } else {
                K = 0
            }
        }
        
        # Set the different slots
        if (!missing(X)) {
            .Object@X = X
        } else {
            .Object@X = matrix(nrow=n, ncol=0)
        }
        if (!missing(V)) {
            .Object@V = V
        } else {
            .Object@V = matrix(nrow=J, ncol=0)
        }
        if (!missing(O_mu)) {
            .Object@O_mu = O_mu
        } else {
            .Object@O_mu = matrix(0, nrow=n, ncol=J)
        }
        if (!missing(O_pi)) {
            .Object@O_pi = O_pi
        } else {
            .Object@O_pi = matrix(0, nrow=n, ncol=J)
        }
        if (!missing(which_X_mu)) {
            .Object@which_X_mu = which_X_mu
        } else if (NCOL(.Object@X)>0) {
            .Object@which_X_mu = seq(NCOL(.Object@X))
        } else {
            .Object@which_X_mu = integer(0)
        }
        if (!missing(which_X_pi)) {
            .Object@which_X_pi = which_X_pi
        } else if (NCOL(.Object@X)>0) {
            .Object@which_X_pi = seq(NCOL(.Object@X))
        } else {
            .Object@which_X_pi = integer(0)
        }
        if (!missing(which_V_mu)) {
            .Object@which_V_mu = which_V_mu
        } else if (NCOL(.Object@V)>0) {
            .Object@which_V_mu = seq(NCOL(.Object@V))
        } else {
            .Object@which_V_mu = integer(0)
        }
        if (!missing(which_V_pi)) {
            .Object@which_V_pi = which_V_pi
        } else if (NCOL(.Object@V)>0) {
            .Object@which_V_pi = seq(NCOL(.Object@V))
        } else {
            .Object@which_V_pi = integer(0)
        }
        if (!missing(W)) {
            .Object@W = W
        } else {
            .Object@W = matrix(0, nrow=n , ncol=K)
        }
        if (!missing(beta_mu)) {
            .Object@beta_mu = beta_mu
        } else {
            .Object@beta_mu = matrix(0, nrow=length(.Object@which_X_mu) , ncol=J)
        }
        if (!missing(beta_pi)) {
            .Object@beta_pi = beta_pi
        } else {
            .Object@beta_pi = matrix(0, nrow=length(.Object@which_X_pi) , ncol=J)
        }
        if (!missing(gamma_mu)) {
            .Object@gamma_mu = gamma_mu
        } else {
            .Object@gamma_mu = matrix(0, nrow=length(.Object@which_V_mu) , ncol=n)
        }
        if (!missing(gamma_pi)) {
            .Object@gamma_pi = gamma_pi
        } else {
            .Object@gamma_pi = matrix(0, nrow=length(.Object@which_V_pi) , ncol=n)
        }
        if (!missing(alpha_mu)) {
            .Object@alpha_mu = alpha_mu
        } else {
            .Object@alpha_mu = matrix(0, nrow=K , ncol=J)
        }
        if (!missing(alpha_pi)) {
            .Object@alpha_pi = alpha_pi
        } else {
            .Object@alpha_pi = matrix(0, nrow=K , ncol=J)
        }
        if (!missing(phi)) {
            .Object@phi = phi
        } else {
            .Object@phi = rep(1,J)
        }
        validObject(.Object) # call of the inspector
        return(.Object)
    }
)

#' Initialize an object of class zinb_model
#' @export
zinb_model <- function(...) {
    new(Class="zinb_model", ...)
}

setMethod("getN", "zinb_model",
          function(object) {
              return(NROW(object@X))
          }
)

setMethod("getJ", "zinb_model",
          function(object) {
              return(NROW(object@V))
          }
)

setMethod("getK", "zinb_model",
          function(object) {
              return(NCOL(object@W))
          }
)

setMethod("getMu", "zinb_model",
          function(object) {
              return(exp(object@X[,object@which_X_mu] %*% object@beta_mu + t(object@V[,object@which_V_mu] %*% object@gamma_mu) + object@O_mu))
          }
)

setMethod("getPi", "zinb_model",
          function(object) {
              return(binomial()$linkinv(object@X[,object@which_X_mu] %*% object@beta_mu + t(object@V[,object@which_V_mu] %*% object@gamma_mu) + object@O_mu))
          }
)

setMethod("getPhi", "zinb_model",
          function(object) {
              return(object@phi)
          }
)

setMethod(
    f="simulateZINB",
    signature="zinb_model",
    definition=function(object, seed, no_cores) {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
            runif(1)
        if (missing(seed)) 
            RNGstate <- get(".Random.seed", envir = .GlobalEnv)
        else {
            R.seed <- get(".Random.seed", envir = .GlobalEnv)
            set.seed(seed)
            RNGstate <- structure(seed, kind = as.list(RNGkind()))
            on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        }
        if (missing(no_cores)) {
            no_cores=1
        }
        mu = getMu(object)
        pi = getPi(object)
        phi = getPhi(object)
        n = getN(object)
        J = getJ(object)
        
        # Simulate negative binomial with the mean matrix and dispersion
        # parameters
        data.nb <- matrix(unlist( parallel::mclapply(seq(n*J), function(i) { rnbinom(1, mu = mu[i] , size = phi[ceiling(i/n)]) }, mc.cores = no_cores) ), nrow = n )
        
        # Simulate the binary dropout matrix. "1" means that a dropout (zero) is
        # observed instead of the value
        data.dropout <- matrix(unlist( parallel::mclapply(seq(n*J), function(i) { rbinom( 1 , size =1 , prob = pi[i] ) }, mc.cores = no_cores) ), nrow = n )
        
        # Matrix of zero-inflated counts
        counts <- data.nb * (1 - data.dropout)
        
        # Fraction of zeros in the matrix
        zero.fraction <- sum ( counts == 0) / (n*J)    
        
        ret = list ( counts = counts , data.nb = data.nb , data.dropouts = data.dropout , zero.fraction = zero.fraction )
        attr(ret, "seed") <- RNGstate
        ret
    }
)
