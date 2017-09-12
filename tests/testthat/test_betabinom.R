context("betabinom")
test_that("example on beta-binomial", {
  library(emdbook)
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
  theta <- runif(nrow(Y), 1, 100)
  prob <- rnorm(nrow(Y),.5,.01)
  ase.cts <- matrix(rbetabinom(prod(dim(Y)), prob=prob,
                               size=Y, theta=rep(theta,ncol(Y))),
                    nrow=nrow(Y))
  idx <- 1:5
  idx2 <- which(condition == "b")
  theta[idx] <- 100
  ase.cts[idx,idx2] <- matrix(rbetabinom(length(idx)*length(idx2), prob=.75,
                                         size=Y[idx,idx2], theta=100),
                              nrow=length(idx))
  betabinom.log.lik <- function(y, x, beta, param, offset) {
    xbeta <- x %*% beta
    p.hat <- (1+exp(-xbeta))^-1
    dbetabinom(y, prob=p.hat, size=param[-1], theta=param[1], log=TRUE)
  }
  theta.hat.0 <- 100 # rough estimate of dispersion
  param <- cbind(theta.hat.0, Y)
  fit.mle <- apeglm(Y=ase.cts, x=x,
                      log.lik=betabinom.log.lik,
                      param=param,
                      no.shrink=TRUE,
                      log.link=FALSE)
  theta.hat <- bbEstDisp(success=ase.cts, size=Y,
                         x=x, beta=fit.mle$map,
                         minDisp=1, maxDisp=500)
  coef <- 2
  mle <- cbind(fit.mle$map[,coef], fit.mle$se[,coef])
  param <- cbind(theta.hat, Y)
  fit2 <- apeglm(Y=ase.cts, x=x,
                 log.lik=betabinom.log.lik,
                 param=param,
                 coef=coef,
                 mle=mle,
                 log.link=FALSE)
  expect_error(Y=ase.cts, x=x,
               log.lik=betabinom.log.lik,
               param=param,
               coef=1,
               mle=mle,
               log.link=FALSE)
  expect_error(Y=ase.cts, x=x,
               log.lik=betabinom.log.lik,
               param=param,
               coef=2,
               log.link=FALSE)
  expect_error(Y=ase.cts, x=x,
               log.lik=betabinom.log.lik,
               param=param,
               coef=2,
               mle = mle,
               log.link=TRUE)
})
