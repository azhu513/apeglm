context("methods")
test_that("alternative methods give same result", {

  n.per.group <- 5
  n <- n.per.group * 2
  m <- 1000
  condition <- factor(rep(letters[1:2], each=n.per.group))
  x <- model.matrix(~condition)
  beta.sd <- 2
  beta.cond <- rnorm(m, 0, beta.sd)
  beta.intercept <- runif(m, 2, 6)
  beta.mat <- cbind(beta.intercept, beta.cond)
  mu <- exp(t(x %*% t(beta.mat)))
  Y <- matrix(rnbinom(m*n, mu=mu, size=1/.1), ncol=n)
  param <- matrix(0.1, nrow=m, ncol=1)
  offset <- matrix(0, nrow=m, ncol=n)
  weights <- matrix(1, nrow=m, ncol=n)
  weights[,c(1,n.per.group+1)] <- 0.01

  # method="general"
  system.time({
    fit <- apeglm(Y=Y, x=x, log.lik=logLikNB, offset=offset, param=param, coef=2)
  })
  #plot(beta.cond, fit$map[,2])

  # method="nbinomR"
  system.time({
    fit.fast <- apeglm(Y=Y, x=x, log.lik=NULL, offset=offset, param=param, coef=2,
                       method="nbinomR")
  })
  expect_equal(fit$map[,1], fit.fast$map[,1], tolerance=1e-3)
  expect_equal(fit$map[,2], fit.fast$map[,2], tolerance=1e-3)
  expect_true(sum(is.na(fit.fast$sd[,1])) == 0)
  
  #plot(fit$map[,1], fit.fast$map[,1])
  #plot(fit$map[,2], fit.fast$map[,2])
  
  # test no offset specified
  fit.fast <- apeglm(Y=Y, x=x, log.lik=NULL, param=param, coef=2, method="nbinomR")

  # test with weights
  fit.w <- apeglm(Y=Y, x=x, log.lik=logLikNB,
                  weights=weights, offset=offset, param=param, coef=2)
  
  fit.w.fast <- apeglm(Y=Y, x=x, log.lik=NULL,
                       weights=weights, offset=offset, param=param, coef=2,
                       method="nbinomR")
  #plot(fit$map[,2], fit.w$map[,2])
  #plot(fit.w$map[,2], fit.w.fast$map[,2])

})
