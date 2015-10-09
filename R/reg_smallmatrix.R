#' Log-likelihood function of the zero-inflated negative binomial model
#' 
#' This function computes the log-likelihood of a standard regression model
#' 
#' This is a (hopefully) memory-efficient implementation of the log-likelihood of a 
#' zero-inflated negative binomial regression model.
#' In this attempt, the design matrices don't have n*J rows, but n and J, respectively.
#' The computation is a bit slower, but the memory usage should be much smaller for
#' large J and n.
#' 
#' @param parms a vector of parameters: should contain the values of beta, followed by those of alpha, followed by the log(1/phi)
#' @param Y the data matrix (genes in rows, cells in columns)
#' @param X the design matrix for the regression on mu (n x k_X)
#' @param W the design matrix for the regression on pi (J x k_W)
#' @param kx the number of beta parameters
#' @param kw the number of alpha parameters
#' @param offsetx the offset for the regression on X
#' @param offsetw the offset for the regression on W
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())
loglik_small <- function(parms, Y, X, W, kx, kw, offsetx, offsetw, linkobj) {
  
  beta <- parms[1:kx]
  mu <- t(exp(X %*% beta + offsetx))
  
  alpha <- matrix(parms[(kx + 1):(kw+kx)], nrow=ncol(W), ncol=n, byrow=TRUE)
  pi <- linkobj$linkinv(W %*% alpha + offsetw)
  
  # for now just one phi
  theta <- exp(parms[(kw + kx) + 1])
  
  loglik0 <- log(pi + exp(log(1 - pi) + suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE))))
  loglik1 <- log(1 - pi) + suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
  
  return(sum(loglik0[which(Y==0)]) + sum(loglik1[which(Y>0)]))
}
