context("intercept")
test_that("intercept issue is solved", {

  set.seed(1)
  n.per.group <- 4
  n <- n.per.group * 2
  m <- 100
  condition <- factor(rep(letters[1:2], each=n.per.group))
  w <- rnorm(n)
  # with this design, no sample has the intercept coefficient alone
  x <- model.matrix(~w + condition)
  beta.sd <- 2
  beta.cond <- rnorm(m, 0, beta.sd)
  beta.intercept <- runif(m, 2, 6)
  beta.mat <- cbind(beta.intercept, rep(0, m), beta.cond)
  mu <- exp(t(x %*% t(beta.mat)))
  Y <- matrix(rnbinom(m*n, mu=mu, size=1/.1), ncol=n)
  param <- matrix(0.1, nrow=m, ncol=1)
  offset <- matrix(0, nrow=m, ncol=n)
  fit <- apeglm(Y=Y, x=x, log.lik=logLikNB, offset=offset, param=param, coef=2)

})
