#' Solve ridge regression or logistic regression problems
#'
#' This function solves a regression or logistic regression problem regularized
#' by a L2 or weighted L2 penalty. Contrary to \code{lm.ridge} or \code{glmnet},
#' it works for any number of predictors.
#' @param x a matrix of covariates, one sample per row, one covariate per
#'   column.
#' @param y a matrix of responses (continuous for regression, 0/1 binary for
#'   logistic regression), one sample per row, one feature per column.
#' @param P a logical matrix of the same size of \code{y} indicating whether
#'   value should be used for the regression (useful to "mask" NA values).
#' @param beta an initial solution where optimization starts (null vector by
#'   default)
#' @param epsilon a scalar or vector of regularization parameters (default
#'   \code{1e-6})
#' @param family a string to choose the type of regression (default
#'   \code{family="gaussian"})
#' @param offset scalar or matrix of offsets (default 0); if a matrix, it must
#'   have the same size of \code{y}
#' @return A vector solution of the regression problem
#' @details When \code{family="gaussian"}, we solve the ridge regression problem
#'   that finds the \eqn{\beta} that minimizes: \deqn{||y - x \beta||^2 +
#'   \epsilon||\beta||^2/2 .} When \code{family="binomial"} we solve the ridge
#'   logistic regression problem \deqn{min \sum_i [-y_i (x \beta)_i +
#'   log(1+exp(x\beta)_i)) ] + \epsilon||\beta||^2/2 .} When \code{epsilon} is a
#'   vector of size equal to the size of \code{beta}, then the penalty is a
#'   weighted L2 norm \eqn{\sum_i \epsilon_i \beta_i^2 / 2}.
#' @export
solveRidgeRegression <- function(x, y,
                                 P = matrix(TRUE, nrow=NROW(y), ncol=NCOL(y)),
                                 beta=rep(0, NCOL(x) * NROW(y)),
                                 epsilon=1e-6, family=c("gaussian","binomial"),
                                 offset=0) {

    family <-  match.arg(family)

    # loglik
    f <- if (family == "gaussian") {
        function(b, P) {
            bb <- matrix(b, ncol=NROW(y), nrow=NCOL(x))
            eta <- t(x %*% bb) + offset
            yeta <- (eta - y)
            yeta[!P] <- 0
            l <- sum((yeta)^2)/2
            l + sum(epsilon*b^2)/2
        }
    } else if (family == "binomial") {
        function(b, P) {
            bb <- matrix(b, ncol=NROW(y), nrow=NCOL(x))
            eta <- t(x %*% bb) + offset
            yeta <- -y*eta + copula::log1pexp(eta)
            yeta[!P] <- 0
            l <- sum(yeta)
            l + sum(epsilon*b^2)/2
        }
    }

    # gradient of loglik
    g <- if (family == "gaussian") {
        function(b, P) {
            bb <- matrix(b, ncol=NROW(y), nrow=NCOL(x))
            eta <- t(x %*% bb) + offset
            yeta <- (eta - y)
            yeta[!P] <- 0
            l <- t(x) %*% t(yeta)
            as.vector(l + epsilon*bb)
        }
    } else if (family == "binomial") {
        function(b, P) {
            bb <- matrix(b, ncol=NROW(y), nrow=NCOL(x))
            eta <- t(x %*% bb) + offset
            yeta <- (-y + 1/(1+exp(-eta)))
            yeta[!P] <- 0
            l <- t(x) %*% t(yeta)
            l + epsilon*b
        }
    }

    # optimize
    m <- optim(par=beta, fn=f, gr=g, P=P, control=list(trace=0), method="BFGS")
    matrix(m$par, nrow=NCOL(x))

}
