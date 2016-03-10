library(sBIC)
library(flexmix)

context("Testing Binomial Mixtures")

test_that("NO TESTS YET IMPLEMENTED", {
  expect_true(FALSE)
})
#
# rmixedbinom = function(n, alpha, theta, size) {
#   ## This function generates ramdom sample from Binomial mixture.
#   ##
#   ## n:      sample size.
#   ## alpha:  mixing proprotions.
#   ## theta:  mixing parameters.
#   ## size:   size parameter for Binomial mixture.
#
#   m = length(theta)
#   data = c()
#   nindex = as.numeric(rmultinom(1, n, alpha))
#   for (i in 1:m) {
#     data = c(data, rbinom(nindex[i], size, theta[i]))
#   }
#   data
# }
#
# ## set.seed(10)  ## run1
# set.seed(11)
# N = 30
# nsim = 200
#
# k.max = 6 # k.max = max # components to be considered
#
# k = 4
# alpha = rep(1 / k, k)
# alpha
# theta = seq(0.2, 0.8, length = k)
# theta
#
# nn <- c(50, 100, 250, 500, 750, 1000)
# nn <- c(50, 100, 200, 300, 400, 500)
# nn <- c(50, 200, 500)
# print(nn)
#
# bicSingular <-
#   bicSingularRouss <- bicRegular <- matrix(0, length(nn), k.max)
# dimnames(bicSingular) <-
#   dimnames(bicSingularRouss) <- dimnames(bicRegular) <-
#   list(as.character(nn), as.character(1:k.max))
# bicS.res <- bicSRouss.res <- bicR.res <- numeric(nsim)
#
# for (n in nn) {
#   print(paste("n =", n))
#   for (j in 1:nsim) {
#     print(j)
#     Y = rmixedbinom(n, alpha, theta, N)
#     X = cbind(Y, N - Y)
#     m <- sBIC(X, BinomialMixtures(k.max))
#
#     reg.BIC <- m$regBIC
#
#     BIC.val <- m$L
#     ## not on log scale yet
#     bicS.res[j] = which(BIC.val == max(BIC.val))
#
#     bicR.res[j] =  which(reg.BIC == max(reg.BIC))
#
#     m <- sBIC(X, BinomialMixtures(k.max, rousseau = T))
#
#     BIC.val <- m$L
#     ## not on log scale yet
#     bicSRouss.res[j] = which(BIC.val == max(BIC.val))
#   }
#
#   bicS.res = factor(bicS.res, levels = seq(1:k.max))
#   bicSRouss.res = factor(bicSRouss.res, levels = seq(1:k.max))
#   bicR.res = factor(bicR.res, levels = seq(1:k.max))
#
#   bicSingular[as.character(n), ] <- table(bicS.res)
#   bicSingularRouss[as.character(n), ] <- table(bicSRouss.res)
#   bicRegular[as.character(n), ] <- table(bicR.res)
#
#   bicSingular[as.character(n), ] <- table(bicS.res)
#   bicSingularRouss[as.character(n), ] <- table(bicSRouss.res)
#   bicRegular[as.character(n), ] <- table(bicR.res)
#   cat(
#     n,
#     "\t",
#     "BICS",
#     table(bicS.res),
#     "\t",
#     "Rouss",
#     table(bicSRouss.res),
#     "\t",
#     "BIC",
#     table(bicR.res),
#     "\n"
#   )
# }
#
# ## 50 	 BICS 0 1 88 11 0 0 0 	 BIC 0 9 86 5 0 0 0
# ## 100 	 BICS 0 0 70 30 0 0 0 	 BIC 0 0 88 12 0 0 0
# ## 200 	 BICS 0 0 45 55 0 0 0 	 BIC 0 0 74 26 0 0 0
# ## 300 	 BICS 0 0 21 79 0 0 0 	 BIC 0 0 56 44 0 0 0
# ## 400 	 BICS 0 0 12 87 1 0 0 	 BIC 0 0 39 61 0 0 0
