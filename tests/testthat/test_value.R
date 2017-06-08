context("Test numerical correctness of functions.")
set.seed(13124)

BiocParallel::register(BiocParallel::SerialParam())

test_that("Estimates are reasonable when data is Poisson", {
    counts <- matrix(rpois(10000, lambda=50), nrow=100, ncol=100)
    m1 <- zinbFit(counts, commondispersion = TRUE)
    expect_true(all(getPhi(m1) < 1e-4))

    m2 <- zinbFit(counts, commondispersion = FALSE)
    expect_true(all(getPhi(m2) < 1e-4))

    expect_equivalent(round(getMu(m1), 2), round(getMu(m2), 2))
    expect_true(abs(mean(getMu(m1)) - 50) < 1)
    expect_true(mean(getPi(m1)) < 1e-2)
})

test_that("Estimates are reasonable when data is Negative Binomial", {
    counts <- matrix(rnbinom(10000, mu=50, size = 10), nrow=100, ncol=100)

    m1 <- zinbFit(counts, commondispersion = TRUE)

    expect_true(abs(mean(getMu(m1)) - 50) < 1)
    expect_true(abs(mean(getTheta(m1)) - 10) < 1)
    expect_true(mean(getPi(m1)) < 1e-2)
})

test_that("solveRidgeRegression gives approximately expected values", {
    set.seed(744747)
    n <- 2
    J <- 10
    Y <- matrix(rnbinom(1000, mu=5, size = 1), nrow=n, ncol=J)
    V <- matrix(c(rep(1, 10), rnorm(10)), ncol=2)
    m1 <- zinbFit(t(Y), V=V, commondispersion = TRUE)

    P <- Y > 0
    L <- matrix(0, nrow=n, ncol=J)
    L[P] <- log(Y[P])

    ## gamma_mu
    Xbeta_mu <- getX_mu(m1) %*% getBeta_mu(m1)
    V <- getV_mu(m1)
    X <- getX_mu(m1)
    epsilon <- getEpsilon_gamma_mu(m1)

    gamma_mu <- matrix(c(1.553248, 0.1422652, 1.80530920,
                         -0.04269995), ncol=2)

    y <- matrix(0, nrow=n, ncol=J)
    y[P] <- L[P] - Xbeta_mu[P]

    gm <- solveRidgeRegression(x=V,
                             y=y,
                             P=P,
                             epsilon = epsilon,
                             family="gaussian")

    expect_equal(round(gamma_mu, 1), round(gm, 1))

    ## beta_mu

    tVgamma_mu <- t(getV_mu(m1) %*% gamma_mu)
    y <- matrix(0, nrow=n, ncol=J)
    y[P] = L[P] - tVgamma_mu[P]

    bm <- solveRidgeRegression(x = getX_mu(m1),
                               y = t(y),
                               P = t(P),
                               epsilon = getEpsilon_beta_mu(m1),
                               family="gaussian")

    beta_mu <- matrix(c(-0.3187399, 0.3596111, 1.091645, 1.490594, 0.02116552,
                        -0.2157922, -0.4769161, 0.5874102, -0.2920499, -1.036666),
                      ncol=10)

    expect_equal(round(beta_mu, 1), round(bm, 1))

    ## check zinbInitialize
    m0 <- zinbModel(n=NROW(Y), J=NCOL(Y), X=getX_mu(m1), V=getV_mu(m1))
    m2 <- zinbInitialize(m0, Y)
    bm <- getBeta_mu(m2)
    bp <- getBeta_pi(m2)
    gm <- getGamma_mu(m2)
    gp <- getGamma_pi(m2)

    beta_mu <- matrix(c(-0.423503, 0.3018034, 0.8690388, 1.345368, -0.02791817,
                        -0.195978, -0.5277613, 0.4950754, -0.4748292, -1.210463),
                      ncol=10)
    beta_pi <- matrix(c(-3.854019, 5.379042, -3.841036, -3.84956, -3.860157,
                        -3.867754, -3.859962, -3.855389, -3.845423, -3.846413),
                      ncol=10)
    gamma_mu <- matrix(c(1.6731321, 0.1998334, 1.91071538, 0.02021601), ncol=2)
    gamma_pi <- matrix(c(-1.702345382, -0.008373722, -9.315069837, 0.001405084),
                       ncol=2)

    expect_equal(round(beta_mu, 1), round(bm, 1))
    expect_equal(round(gamma_mu, 1), round(gm, 1))
    expect_equal(round(beta_pi), round(bp))
    expect_equal(round(gamma_pi, 1), round(gp, 1))

})
