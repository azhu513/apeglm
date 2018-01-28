context("methods")
test_that("alternative methods give same result", {

  set.seed(1)
  n.per.group <- 10
  n <- n.per.group * 2
  m <- 1000
  condition <- factor(rep(1:2, each=n.per.group))
  batch <- factor(rep(rep(1:2, each=n.per.group/2),2))
  x <- model.matrix(~batch + condition)
  beta.cond <- rnorm(m, 0, 2)
  beta.batch <- rnorm(m, 0, 1)
  beta.intercept <- runif(m, 2, 6)
  beta.mat <- cbind(beta.intercept, beta.batch, beta.cond)
  q <- exp(t(x %*% t(beta.mat)))
  sf <- 2^rep(seq(-.5,.5,length=n.per.group),2)
  mu <- t(t(q) * sf)
  Y <- matrix(rnbinom(m*n, mu=mu, size=1/.1), ncol=n)
  param <- matrix(0.1, nrow=m, ncol=1)

  offset <- matrix(log(sf), nrow=m, ncol=n, byrow=TRUE)
  weights <- matrix(1, nrow=m, ncol=n)
  weights[,c(1,n.per.group+1)] <- 0.01

  # method="general"
  coef <- 3
  system.time({
    fit <- apeglm(Y=Y, x=x, log.lik=logLikNB, offset=offset, param=param, coef=coef)
  })

  # method="nbinomR"
  system.time({
    fitR <- apeglm(Y=Y, x=x, log.lik=NULL, offset=offset, param=param, coef=coef,
                       method="nbinomR")
  })
  expect_equal(fit$map[,coef], fitR$map[,coef], tolerance=1e-3)
  expect_equal(fit$sd[,coef], fitR$sd[,coef], tolerance=1e-3)
  
  # method="nbinomR" with test no offset specified
  fitNoOff <- apeglm(Y=Y, x=x, log.lik=NULL, param=param, coef=coef,
                     method="nbinomR")

  # pretty fast in C++, only returns MAP coefficients
  system.time({
    fitC <- apeglm(Y=Y, x=x, log.lik=NULL, offset=offset, param=param, coef=coef,
                   method="nbinomC")
  })
  expect_equal(fit$map[,coef], fitC$map[,coef], tolerance=1e-3)
  #expect_equal(fitC$diag[,"value"], fitC$diag[,"valueR"], tolerance=1e-3)

  # C++ to fit the MAP, then estimate posterior SD in R
  system.time({
    fitCR <- apeglm(Y=Y, x=x, log.lik=NULL, offset=offset, param=param, coef=coef,
                    method="nbinomCR")
  })
  expect_equal(fit$sd[,coef], fitCR$sd[,coef], tolerance=1e-3)

  # test with weights
  fit.w <- apeglm(Y=Y, x=x, log.lik=logLikNB,
                  weights=weights, offset=offset, param=param, coef=coef)  
  fitR.w <- apeglm(Y=Y, x=x, log.lik=NULL,
                   weights=weights, offset=offset, param=param, coef=coef,
                   method="nbinomR")
  fitCR.w <- apeglm(Y=Y, x=x, log.lik=NULL,
                    weights=weights, offset=offset, param=param, coef=coef,
                    method="nbinomCR")
  expect_equal(fit.w$map[,coef], fitR.w$map[,coef], tolerance=1e-3)
  expect_equal(fit.w$sd[,coef], fitR.w$sd[,coef], tolerance=1e-3)
  expect_equal(fit.w$map[,coef], fitCR.w$map[,coef], tolerance=1e-3)
  expect_equal(fit.w$sd[,coef], fitCR.w$sd[,coef], tolerance=1e-3)
  
})
