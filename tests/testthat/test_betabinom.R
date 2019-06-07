context("betabinom")
test_that("example on beta-binomial", {

  library(emdbook)
  
  set.seed(1)
  n.per.group <- 20
  n <- 2 * n.per.group
  condition <- factor(rep(1:2,each=n.per.group))
  m <- 100
  cts <- matrix(rnbinom(m*n, mu=200, size=1/.1), ncol=n)
  theta <- runif(m,1,100)
  
  x <- model.matrix(~condition)
  beta.sd <- 2
  beta.cond <- rnorm(m, 0, beta.sd)
  beta.intercept <- rnorm(m, 0, .1)
  beta.mat <- cbind(beta.intercept, beta.cond)
  prob <- (1+exp(-(beta.mat %*% t(x))))^-1
  ase.cts <- matrix(rbetabinom(prod(dim(cts)), prob=prob,
                               size=cts, theta=rep(theta,ncol(cts))),
                    nrow=nrow(cts))
  
  betabinom.log.lik <- function(y, x, beta, param, offset) {
    xbeta <- x %*% beta
    p.hat <- (1+exp(-xbeta))^-1
    dbetabinom(y, prob=p.hat, size=param[-1], theta=param[1], log=TRUE)
  }
  
  theta.hat.0 <- 100 # rough estimate of dispersion
  param <- cbind(theta.hat.0, cts)
  fit.mle <- apeglm(Y=ase.cts, x=x, log.lik=betabinom.log.lik,
                    param=param, no.shrink=TRUE, log.link=FALSE)
  
  theta.hat <- bbEstDisp(success=ase.cts, size=cts,
                         x=x, beta=fit.mle$map,
                         minDisp=1, maxDisp=500)
  
  # we can also ask for standard errors for log of dispersion
  theta.fit <- bbEstDisp(success=ase.cts, size=cts,
                         x=x, beta=fit.mle$map,
                         minDisp=1, maxDisp=500, se=TRUE)

  #plot(beta.cond, fit.mle$map[,2])
  #plot(theta, theta.hat)
  
  coef <- 2
  mle <- cbind(fit.mle$map[,coef], fit.mle$sd[,coef])
  param <- cbind(theta.hat, cts)

  fit <- apeglm(Y=ase.cts, x=x, log.lik=betabinom.log.lik,
                param=param, coef=coef, mle=mle, log.link=FALSE)

  #plot(beta.cond, fit$map[,2])
  system.time({
  fitC <- apeglm(Y=ase.cts, x=x, method="betabinC",
                 param=param, coef=coef, mle=mle)
  })
  
  fitCstar <- apeglm(Y=ase.cts, x=x, method="betabinC*",
                     param=param, coef=coef, mle=mle)

  fitR <- apeglm(Y=ase.cts, x=x, method="betabinR",
                 param=param, coef=coef, mle=mle)

  fitCR <- apeglm(Y=ase.cts, x=x, method="betabinCR",
                  param=param, coef=coef, mle=mle)
  
  expect_equal(fit$map[,2], fitC$map[,2], tolerance=1e-3)
  expect_equal(fit$map[,2], fitCstar$map[,2], tolerance=1e-3)
  expect_equal(fit$map[,2], fitR$map[,2], tolerance=1e-3)
  expect_equal(fit$map[,2], fitCR$map[,2], tolerance=1e-3)

  expect_equal(fit$sd[,2], fitR$sd[,2], tolerance=1e-3)
  expect_equal(fit$sd[,2], fitCR$sd[,2], tolerance=1e-3)
  
})
