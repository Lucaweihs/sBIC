library(sBIC)

context("Testing ReducedRankRegressions")

test_that("NOT ENOUGH TESTS IMPLEMENTED", {
  expect_true(FALSE)
})

# rortho = function(n) {
#   return(qr.Q(qr(array(runif(n), dim = c(n, n)))))
# }
#
# set.seed(1234)
#
# for (ii in 1:10){
#   print(ii)
#   sims = 500
#
#   M = 15
#   N = 10
#   nn = c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500)
#   n = nn[ii]
#   r = 4
#   sing.vals = (r:1) / 4 + 1 / 4
#
#   rank.bic = factor(rep(0, sims), levels = 0:min(M, N))
#   rank.bics = factor(rep(0, sims), levels = 0:min(M, N))
#   ## +1 for 0:min(M,N)
#
#   for (i in 1:sims) {
#     B = rortho(N)[,1:r] %*% diag(sing.vals) %*% rortho(M)[1:r,]
#     X = matrix(rnorm(M * n), M, n)
#     Y = 1 * B %*% X + matrix(rnorm(N * n), N, n)
#
#     maxRank = min(M, N)
#     rrRegressions = ReducedRankRegressions(N, M, maxRank)
#     results = sBIC(list(X = X, Y = Y), rrRegressions)
#
#     rank.bic[i] = which.max(results$regBIC) - 1
#     rank.bics[i] = which.max(results$logL) - 1
#   }
#
#   tt = rbind(table(rank.bic), table(rank.bics))
#   dimnames(tt)[[1]] = c("BIC", "sBIC")
#   print(n)
#   print(tt)
# }