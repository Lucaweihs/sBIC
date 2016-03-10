library(sBIC)

context("Testing Factor Analyses")

test_that("NOT ENOUGH TESTS IMPLEMENTED", {
  expect_true(FALSE)
})

test_that("Learning coefficients are correct", {
  facModels = FactorAnalyses(6, 3)
  trueLearningCoefs = matrix(c(3,0,0,0,
                               9/2,6,0,0,
                               6,29/4,17/2,0,
                               15/2,17/2,19/2,21/2), ncol = 4, byrow = T)
  for (i in 1:4) {
    for (j in 1:i) {
      expect_equal(facModels$learnCoef(i, j)$lambda, trueLearningCoefs[i, j])
      expect_equal(facModels$learnCoef(i, j)$m, 1)
    }
  }
})

# X = read.table("lopeswest.csv", header = TRUE)[, -(1:2)]
# m = ncol(X)
# kmax = 3
# Y = scale(X[-1, ] - X[-n, ])
# facModels = FactorAnalyses(m, kmax)
# sBIC(Y, facModels)

