context("nbinom")
test_that("nbinom cases works", {
  n.per.group <- 5
  n <- n.per.group * 2
  m <- 100
  condition <- factor(rep(letters[1:2], each = n.per.group))
  x <- model.matrix(~condition)
  beta.sd <- 2
  beta.cond <- rnorm(m, 0, beta.sd)
  beta.intercept <- runif(m, 2, 6)
  beta.mat <- cbind(beta.intercept, beta.cond)

  mu <- exp(t(x %*% t(beta.mat)))
  Y <- matrix(rnbinom(m*n, mu=mu, size=1/.1), ncol = n)

  param <- matrix(0.1, nrow = m, ncol = 1)
  offset <- matrix(0, nrow = m, ncol = n)

  res <- apeglm(Y = Y, x = x,
                log.lik = logLikNB,
                param = param,
                offset = offset,
                coef = 2)
})