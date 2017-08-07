#' Simple line-search estimator for dispersion of a beta binomial
#'
#' Uses R's \code{optimize} function to find the maximum likelihood
#' estimate of dispersion for a beta binomial distribution
#' (\code{theta} for the \code{dbetabinom} function in the
#' emdbook package). The counts, size, and beta are matrices,
#' such that each row could be treated as a beta-binomial GLM
#' problem.
#'
#' @param success the observed successes (a matrix)
#' @param size the total trials (a matrix)
#' @param x the design matrix, as many rows as columns of \code{success} and \code{size}
#' @param beta a matrix of MLE coefficients, as many rows as \code{success} and \code{size}
#' @param minDisp the minimum dispersion value
#' @param maxDisp the maximum dispersion value
#'
#' @return a vector of estimated dispersions
#'
#' @examples
#'
#' library(emdbook)
#' n <- 100
#' m <- 100
#' size <- matrix(rnbinom(n*m, mu=100, size=10),ncol=m)
#' success <- matrix(rbetabinom(n*m, prob=.5, size=size, theta=100),ncol=m)
#' x <- matrix(rep(1,m),ncol=1)
#' beta <- matrix(rep(0,n),ncol=1)
#' theta <- bbEstDisp(success, size, x, beta, 1, 500)
#' summary(theta)
#'
#' @importFrom stats optimize
#' @importFrom emdbook dbetabinom
#' 
#' @export
bbEstDisp <- function(success, size, x, beta, minDisp, maxDisp) {
  stopifnot(ncol(success) == nrow(x))
  stopifnot(nrow(success) == nrow(beta))
  stopifnot(ncol(x) == ncol(beta))
  xbeta <- t(x %*% t(beta))
  p.hat <- (1+exp(-xbeta))^-1
  theta.hat <- numeric(nrow(success))
  minld <- log(minDisp)
  maxld <- log(maxDisp)
  for (i in seq_len(nrow(success))) {
    f <- function(logtheta,i) sum(dbetabinom(success[i,], p.hat[i,], size=size[i,],
                                             theta=exp(logtheta), log=TRUE))
    o <- optimize(f, interval=c(minld,maxld), i=i, maximum=TRUE)
    theta.hat[i] <- exp(o$maximum)
  }
  theta.hat
}
