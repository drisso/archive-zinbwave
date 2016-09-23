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

