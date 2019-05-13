# negative log posterior
nbinomFn <- function(beta, x, y, size, weights, offset, sigma, S, no.shrink, shrink, cnst) {
  xbeta <- x %*% beta
  prior <- sum(-beta[no.shrink]^2/(2*sigma^2)) + sum(-log1p(beta[shrink]^2/S^2))
  -sum(weights * (y * xbeta - (y + size) * log(size + exp(xbeta + offset)))) - prior + cnst
}

# gradient of negative log posterior
nbinomGr <- function(beta, x, y, size, weights, offset, sigma, S, no.shrink, shrink, cnst) {
  xbeta <- x %*% beta
  exbetaoff <- exp(xbeta + offset)
  prior <- numeric(length(beta))
  prior[no.shrink] <- -beta[no.shrink]/sigma^2
  prior[shrink] <- -2*beta[shrink]/(S^2 + beta[shrink]^2)
  -t(x) %*% (weights * (y - (y + size) * exbetaoff / (size + exbetaoff))) - prior
}

optimNbinom <- function(init, y, x, param, weights, offset, prior.control,
                        bounds, optim.method) {
  if (is.null(weights)) {
    weights <- 1
  }
  if (is.null(offset)) {
    offset <- 0
  }
  size <- 1/param
  no.shrink <- prior.control$no.shrink
  shrink <- setdiff(seq_along(init), no.shrink)
  sigma <- prior.control$prior.no.shrink.scale
  S <- prior.control$prior.scale
  cnst <- -nbinomFn(init, x, y, size, weights, offset, sigma, S, no.shrink, shrink, 0) - 1
  o <- optim(par=init, fn=nbinomFn, gr=nbinomGr,
             x=x, y=y, size=size,
             weights=weights, offset=offset,
             sigma=sigma, S=S, no.shrink=no.shrink,
             shrink=shrink, cnst=cnst,
             lower=bounds[1], upper=bounds[2],
             hessian=TRUE, method=optim.method)
  o$hessian <- -1 * o$hessian
  o
}

optimNbinomHess <- function(init, y, x, param, weights, offset, prior.control,
                            bounds, optim.method, prefit.conv) {
  if (is.null(weights)) {
    weights <- 1
  }
  if (is.null(offset)) {
    offset <- 0
  }
  size <- 1/param
  no.shrink <- prior.control$no.shrink
  shrink <- setdiff(seq_along(init), no.shrink)
  sigma <- prior.control$prior.no.shrink.scale
  S <- prior.control$prior.scale
  cnst <- -nbinomFn(init, x, y, size, weights, offset, sigma, S, no.shrink, shrink, 0) - 1
  o <- list()
  o$par <- init
  if (prefit.conv == 0) {
    o$hessian <- -1 * optimHess(par=init, fn=nbinomFn, gr=nbinomGr,
                                x=x, y=y, size=size,
                                weights=weights, offset=offset,
                                sigma=sigma, S=S, no.shrink=no.shrink,
                                shrink=shrink, cnst=cnst)
    var.est <- diag(-solve(o$hessian))
  }
  # rows that did not converge in C++ will get another try below
  o$convergence <- 0
  o$counts <- NA
  # if the C++ fit did not converge, or we have negative variance estimates
  # run the optimization in R
  if (prefit.conv != 0 || any(var.est <= 0)) {
    o <- optimNbinom(init=init, y=y, x=x, param=param,
                     weights=weights, offset=offset,
                     prior.control=prior.control,
                     bounds=bounds, optim.method=optim.method)
  }
  # don't want to mix C++ and R log posterior values, which have diff scale
  o$value <- NA
  o
}


######################################################
## the following is for beta-binomial               ##
## written by Josh Zitovsky in Spring semseter 2019 ##
######################################################

# negative log posterior (BB)
betabinFn <- function(beta, x, y, size, theta, weights, sigma, S, no.shrink, shrink, cnst) {
  # prior log-likelihood for shrunk and non-shrunk estimates
  prior <- sum(-beta[no.shrink]^2/(2*sigma^2)) + sum(-log1p(beta[shrink]^2/S^2))
  xbeta <- x %*% beta
  # estimated probabilities for a GLM with a logit link function
  exp_neg_xbeta <- exp(-xbeta)
  p.hat <- (1 + exp_neg_xbeta)^-1
  ptheta <- p.hat * theta
  # vector of likelihoods for each observation (sample)
  observedParts <- suppressWarnings(lgamma(size - y + theta - ptheta) +
                                    lgamma(y + ptheta) -
                                    lgamma(theta - ptheta) -
                                    lgamma(ptheta))
  #joint observed data likelihood
  observed <- sum(observedParts*weights)
  # -L+L0+1, where L is the posterior likelihood and L0 is the
  # posterior lilelihood of the initial or prefit beta values.
  # The addition of the constant normally prevents the returning
  # value from being too large or small.
  return(-prior - observed + cnst)
}

# gradient of negative log posterior (BB)
betabinGr <- function(beta, x, y, size, theta, weights, sigma, S, no.shrink, shrink, cnst) {
  prior <- numeric(length(beta))
  # partial derivatives of the log prior for non-shrunk estimates
  prior[no.shrink] <- -beta[no.shrink]/sigma^2
  # partial derivatives of the log prior for shrunk-estimates
  prior[shrink] <- -2*beta[shrink]/(S^2 + beta[shrink]^2)

  xbeta <- x %*% beta
  exp_neg_xbeta <- exp(-xbeta)
  # estimated probabilities
  p.hat <- (1 + exp_neg_xbeta)^-1
  ptheta <- p.hat * theta
  e1 <- (theta * exp_neg_xbeta)/((1 + exp_neg_xbeta)^2)
  e2 <- suppressWarnings(-digamma(size - y + theta - ptheta) +
                         digamma(y + ptheta) +
                         digamma(theta - ptheta) -
                         digamma(ptheta))
  e3 <- e1 * e2 * weights
  # gradient vector of the data log-likelihood
  observed <- (t(x) %*% e3)
  # gradient vector of the negative log posterior
  return(-observed - prior)
}

  # similar to optimNbinom
  optimBetabin <- function(init, y, x, size, theta, weights, prior.control,
                          bounds, optim.method, cnst=0) {
    if (is.null(weights)) {
      weights <- 1
    }
    no.shrink <- prior.control$no.shrink
    shrink <- setdiff(seq_along(init), no.shrink)
    sigma <- prior.control$prior.no.shrink.scale
    S <- prior.control$prior.scale
    cnst <- -betabinFn(init, x, y, size, theta, weights, sigma, S, no.shrink, shrink, 0) - 1
    o <- optim(par=init, fn=betabinFn, gr=betabinGr,
               x=x, y=y, size=size,
               theta=theta, weights=weights,
               sigma=sigma, S=S, no.shrink=no.shrink,
               shrink=shrink, cnst=cnst,
               lower=bounds[1], upper=bounds[2],
               hessian=TRUE, method=optim.method)
    o$hessian <- -1 * o$hessian
    o
  }

  # similar to optimNbinomHess
  optimBetabinHess <- function(init, y, x, param, weights, prior.control,
                               bounds, optim.method, prefit.conv) {
  if (is.null(weights)) {
    weights <- 1
  }
  # overdispersion
  theta <- param[1]
  # vector of sizes (total counts)
  size <- param[-1]
  # variables to not shrink
  no.shrink <- prior.control$no.shrink
  # variables to shrink
  shrink <- setdiff(seq_along(init), no.shrink)
  sigma <- prior.control$prior.no.shrink.scale
  S <- prior.control$prior.scale
  o <- list()
  o$par <- init
  cnst <- -betabinFn(init, x, y, size, theta, weights, sigma, S, no.shrink, shrink, 0) - 1
  if (prefit.conv == 0) {
    o$hessian <- -1 * optimHess(par=init, fn=betabinFn, gr=betabinGr,
                                x=x, y=y, size=size, theta=theta,
                                weights=weights,
                                sigma=sigma, S=S, no.shrink=no.shrink,
                                shrink=shrink, cnst=cnst)
    var.est <- diag(-solve(o$hessian))
  }
  o$convergence <- 0
  o$counts <- NA

  # if the C++ fit did not converge, or we have
  # negative variance estimates, run the optimization in R
  if (prefit.conv != 0 || any(var.est <= 0)) {
    o <- optimBetabin(init=init, y=y, x=x, size=size,
                     theta=theta, weights=weights,
                     prior.control=prior.control,
                     bounds=bounds, optim.method=optim.method, cnst=cnst)
  }

  # return R value if betabinCR is selected
  # o$value <- optimBetabin(init=init, y=y, x=x, size=size,
  #                         theta=theta, weights=weights,
  #                         prior.control=prior.control,
  #                         bounds=bounds, optim.method=optim.method, adjust=TRUE)$value
  o$value <- NA
  return(o)
}
