# negative log posterior
nbinomFn <- function(beta, x, y, size, weights, offset, sigma, S, no.shrink, shrink, const) {
  xbeta <- x %*% beta
  prior <- sum(-beta[no.shrink]^2/(2*sigma^2)) + sum(-log1p(beta[shrink]^2/S^2))
  -sum(weights * (y * xbeta - (y + size) * log(size + exp(xbeta + offset)))) - prior + const
}

# gradient of negative log posterior
nbinomGr <- function(beta, x, y, size, weights, offset, sigma, S, no.shrink, shrink, const) {
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
  const <- -nbinomFn(init, x, y, size, weights, offset, sigma, S, no.shrink, shrink, 0) - 1
  o <- optim(par=init, fn=nbinomFn, gr=nbinomGr,
             x=x, y=y, size=size,
             weights=weights, offset=offset,
             sigma=sigma, S=S, no.shrink=no.shrink,
             shrink=shrink, const=const,
             lower=bounds[1], upper=bounds[2],
             hessian=TRUE, method=optim.method)
  o$hessian <- -1 * o$hessian
  o
}

optimNbinomHess <- function(init, y, x, param, weights, offset, prior.control,
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
  const <- -nbinomFn(init, x, y, size, weights, offset, sigma, S, no.shrink, shrink, 0) - 1
  -1 * optimHess(par=init, fn=nbinomFn, gr=nbinomGr,
                 x=x, y=y, size=size,
                 weights=weights, offset=offset,
                 sigma=sigma, S=S, no.shrink=no.shrink,
                 shrink=shrink, const=const)
}
