context("methods")
test_that("alternative methods give same result", {

  set.seed(1)
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
  sf <- 2^rep(seq(-.5,.5,length=n.per.group),2)
  offset <- matrix(log(sf), nrow=m, ncol=n, byrow=TRUE)
  weights <- matrix(1, nrow=m, ncol=n)
  weights[,c(1,n.per.group+1)] <- 0.01

  # method="general"
  system.time({
    fit <- apeglm(Y=Y, x=x, log.lik=logLikNB, offset=offset, param=param, coef=2)
  })

  # method="nbinomR"
  system.time({
    fitR <- apeglm(Y=Y, x=x, log.lik=NULL, offset=offset, param=param, coef=2,
                       method="nbinomR")
  })
  expect_equal(fit$map[,1], fitR$map[,1], tolerance=1e-3)
  expect_equal(fit$map[,2], fitR$map[,2], tolerance=1e-3)
  expect_true(sum(is.na(fitR$sd[,1])) == 0)
  
  # test no offset specified
  fitR <- apeglm(Y=Y, x=x, log.lik=NULL, param=param, coef=2, method="nbinomR")

  # test with weights
  fit.w <- apeglm(Y=Y, x=x, log.lik=logLikNB,
                  weights=weights, offset=offset, param=param, coef=2)  
  fitR.w <- apeglm(Y=Y, x=x, log.lik=NULL,
                   weights=weights, offset=offset, param=param, coef=2,
                   method="nbinomR")
  expect_equal(fit.w$map[,2], fitR.w$map[,2], tolerance=1e-3)
  expect_true(sum(is.na(fitR.w$sd[,1])) == 0)

  # pretty fast in C++, only returns MAP coefficients
  system.time({
    fitC <- apeglm(Y=Y, x=x, log.lik=NULL, offset=offset, param=param, coef=2,
                   method="nbinomC")
  })

  expect_equal(fit$map[,1], fitC$map[,1], tolerance=1e-3)
  expect_equal(fit$map[,2], fitC$map[,2], tolerance=1e-3)
  #expect_equal(fitC$diag[,"value"], fitC$diag[,"valueR"], tolerance=1e-3)

  system.time({
    fitCR <- apeglm(Y=Y, x=x, log.lik=NULL, offset=offset, param=param, coef=2,
                       method="nbinomCR")
  })
  expect_equal(fit$sd[,1], fitCR$sd[,1], tolerance=1e-3)
  expect_equal(fit$sd[,2], fitCR$sd[,2], tolerance=1e-3)
  
})
