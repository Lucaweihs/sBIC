library(sBIC)
library(poLCA)

context("Testing LCAs")

test_that("NOT ENOUGH TESTS IMPLEMENTED", {
  expect_true(FALSE)
})

# set.seed(10)
# nsim = 100
# k.max = 6
# nn <- c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 750, 1000)
# phi = c((1:20) / 2, 1e6)
# no.class = 4
# no.var = 8
# no.levels = 2
# p1 <- 0.85
# p2 <- 0.1
# p3 <- 0.2
# prob.mx <- matrix(c(
#   p1,p2,p2,p2,
#   p1,p3,p3,p3,
#   p2,p1,p2,p2,
#   p3,p1,p3,p3,
#   p2,p2,p1,p2,
#   p3,p3,p1,p3,
#   p2,p2,p2,p1,
#   p3,p3,p3,p1), no.var, no.class, byrow=TRUE)
#
# probs = vector("list", no.var)
# for (i in 1:no.var) {
#   probs[[i]] <- matrix(0, no.class, no.levels)
#   probs[[i]][, 1] <- prob.mx[i, ]
#   probs[[i]][, 2] <- 1 - prob.mx[i, ]
# }
#
# P = rep(1/no.class,no.class)
# P = c(5, 10, 15, 75) / 100  ## unequal
#
# bicRegular <- matrix(0,length(nn),k.max)
# dimnames(bicRegular) <- list(as.character(nn),as.character(1:k.max))
# bicR.res <- numeric(nsim)
#
# bicSingular <- vector("list",length(phi))
# bicS.res <- vector("list",length(phi))
# for(j in 1:length(phi)){
#   bicSingular[[j]] <- matrix(0,length(nn),k.max)
#   dimnames(bicSingular[[j]]) <- list(as.character(nn),as.character(1:k.max))
#   bicS.res[[j]] <- numeric(nsim)
# }
#
# for (N in nn) {
#   for (j in 1:nsim) {
#     print(j)
#     X = poLCA.simdata(N = N, probs = probs, P = P)$dat
#
#     m = sBIC(X, LCAs(k.max, no.var, no.levels, phi[1]))
#
#     BIC.val <- m$L
#     reg.BIC <- m$regBIC
#     bicS.res[[1]][j] = which(BIC.val == max(BIC.val))
#     bicR.res[j] =  which(reg.BIC == max(reg.BIC))
#
#     if (length(phi) > 1) {
#       for (jj in 2:length(phi)) {
#         m = sBIC(X, LCAs(k.max, no.var, no.levels, phi[jj]))
#         BIC.val <- m$L
#         bicS.res[[jj]][j] = which(BIC.val == max(BIC.val))
#
#       }
#     }
#   }
#
#   bicR.res = factor(bicR.res, levels = seq(1:k.max))
#   bicRegular[as.character(N), ] <- table(bicR.res)
#
#   for (jj in 1:length(phi)) {
#     bicS.res[[jj]] = factor(bicS.res[[jj]], levels = seq(1:k.max))
#     bicSingular[[jj]][as.character(N), ] <- table(bicS.res[[jj]])
#   }
#
#   cat(N, "\t", "BIC", table(bicR.res), "\t")
#   for (jj in 1:length(phi)) {
#     cat(paste("BICS", phi[jj], sep = ""), table(bicS.res[[jj]]), "\t")
#   }
#   cat("\n")
# }