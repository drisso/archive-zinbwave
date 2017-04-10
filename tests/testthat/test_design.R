context("Test zinbModel class and methods.")
set.seed(13124)

test_that("getX et al work with/without intercept", {

    bio <- gl(2, 3)
    gc <- rnorm(10)
    m <- zinbFit(matrix(10, 10, 6), X=model.matrix(~bio), V=model.matrix(~gc),
                 which_X_pi=1L, which_V_mu=1L)
    
    expect_equal(NCOL(getV_mu(m)), 1)
    expect_equal(NCOL(getV_mu(m, intercept=FALSE)), 0)
    expect_equal(NCOL(getV_pi(m)), 2)
    expect_equal(NCOL(getV_pi(m, intercept=FALSE)), 1)
    expect_equal(NCOL(getX_mu(m)), 2)
    expect_equal(NCOL(getX_mu(m, intercept=FALSE)), 1)
    expect_equal(NCOL(getX_pi(m)), 1)
    expect_equal(NCOL(getX_pi(m, intercept=FALSE)), 0)
    
    
    m <- zinbFit(matrix(10, 10, 6), X=model.matrix(~bio), V=model.matrix(~gc),
                 which_X_pi=2L, which_V_mu=2L)
    
    expect_equal(getV_mu(m, intercept=TRUE), getV_mu(m, intercept=TRUE))
    expect_equal(getX_pi(m, intercept=TRUE), getX_pi(m, intercept=TRUE))
})

test_that("zinbFit works with genewise dispersion", {
    bio <- gl(2, 3)
    counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
    m <- zinbFit(counts, X=model.matrix(~bio), commondispersion = TRUE)
    m <- zinbFit(counts, X=model.matrix(~bio), commondispersion = FALSE)
})

test_that("zinbSim works", {
    a <- zinbModel(n=5, J=10)
    zinbSim(a)
})

test_that("getMu and getPi have the right dimensions", {
    bio <- gl(2, 3)
    counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
    m <- zinbFit(counts, X=model.matrix(~bio), commondispersion = TRUE)

    expect_equal(dim(getMu(m)), c(nSamples(m), nFeatures(m)))
    expect_equal(dim(getLogMu(m)), c(nSamples(m), nFeatures(m)))
    expect_equal(dim(getPi(m)), c(nSamples(m), nFeatures(m)))
    expect_equal(dim(getLogitPi(m)), c(nSamples(m), nFeatures(m)))
    expect_equal(dim(getW(m)), c(nSamples(m), nFactors(m)))
    expect_equal(length(getPhi(m)), nFeatures(m))
    expect_equal(length(getTheta(m)), nFeatures(m))
    expect_equal(length(getZeta(m)), nFeatures(m))
})

test_that("Initialization works", {
    
    ## no arguments specified
    zinbModel()
    
    ## specify W
    W <- matrix(rnorm(10), ncol=2)
    zinbModel(W = W)
    
    ## add more
})