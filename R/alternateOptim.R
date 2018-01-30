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
