context("nbinom")
test_that("nbinom cases works", {

  set.seed(1)
  n.per.group <- 5
  n <- n.per.group * 2
  m <- 100
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

  # a default run:
  system.time({
    fit <- apeglm(Y=Y, x=x, log.lik=logLikNB, offset=offset, param=param, coef=2)
  })
  #plot(beta.cond, fit$map[,2])
  
  # other interval types:
  fit <- apeglm(Y=Y, x=x,
                log.lik=logLikNB,
                offset=offset,
                param=param,
                coef=2,
                interval.type="HPD")

  fit <- apeglm(Y=Y, x=x,
                log.lik=logLikNB,
                offset=offset,
                param=param,
                coef=2,
                interval.type="credible")
  
  #################
  ## some errors ##
  #################

  # missing 'offset'
  expect_error(apeglm(Y=Y, x=x,
                      log.lik=logLikNB,
                      param=param,
                      coef=2), "offset")
  
  # error for missing param
  expect_error(apeglm(Y=Y, x=x,
                log.lik=logLikNB,
                offset=offset,
                coef=2))

  # error for coef not > 1
  expect_error(apeglm(Y=Y, x=x,
                      log.lik=logLikNB,
                      param=param,
                      offset=offset,
                      coef=1))
  
  # error for threshold not > 0
  expect_error(apeglm(Y=Y, x=x,
                log.lik=logLikNB,
                param=param,
                offset=offset,
                coef=2, threshold=-1))

  # error for threshold & param.sd at same time
  expect_error(apeglm(Y=Y, x=x,
                log.lik=logLikNB,
                param=param,
                offset=offset,
                coef=2, threshold=0, param.sd=0.1))

  # error for wrong shape of offset
  offset <- matrix(0, nrow=n, ncol=m)
  expect_error(apeglm(Y=Y, x=x,
                log.lik=logLikNB,
                param=param,
                offset=offset,
                coef=2)) 
  
})
