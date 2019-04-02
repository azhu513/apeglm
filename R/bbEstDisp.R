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
#' @param weights the weights (1 or a matrix)
#' @param x the design matrix, as many rows as columns of \code{success} and \code{size}
#' @param beta a matrix of MLE coefficients, as many rows as \code{success} and \code{size}
#' @param minDisp the minimum dispersion value
#' @param maxDisp the maximum dispersion value
#' @param se logical, whether to return standard error estimate on the log of
#' the dispersion (theta)
#'
#' @return a vector of estimated dispersions (theta). if \code{se=TRUE} a matrix
#' with columns: the vector of estimated dispersions and the standard
#' errors for the log of the estimated dispersions
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
#' theta <- bbEstDisp(success=success, size=size, x=x, beta=beta, minDisp=1, maxDisp=500)
#' summary(theta)
#'
#' # with standard error estimates on log of dispersion
#' fit <- bbEstDisp(success=success, size=size, x=x, beta=beta, minDisp=1, maxDisp=500, se=TRUE)
#' plot(fit[1:20,"theta"], ylim=c(0,500), ylab="theta-hat")
#' log.theta <- log(fit[1:20,"theta"])
#' log.theta.se <- fit[1:20,"se"]
#' segments(1:20, exp(log.theta - 2 * log.theta.se),
#'          1:20, exp(log.theta + 2 * log.theta.se))
#' abline(h=100,col="red")
#' 
#'
#' @importFrom stats optimize
#' @importFrom emdbook dbetabinom
#' 
#' @export
bbEstDisp <- function(success, size, weights=1, x, beta, minDisp, maxDisp, se=FALSE) {
  stopifnot(ncol(success) == nrow(x))
  stopifnot(nrow(success) == nrow(beta))
  stopifnot(ncol(x) == ncol(beta))
  stopifnot(all(weights >= 0))
  xbeta <- t(x %*% t(beta))
  p.hat <- (1+exp(-xbeta))^-1
  theta.hat <- numeric(nrow(success))
  if (se) {
    se.vec <- numeric(nrow(success))
  }
  minld <- log(minDisp)
  maxld <- log(maxDisp)
  for (i in seq_len(nrow(success))) {
    f <- function(logtheta,i) sum(weights *
                                  dbetabinom(success[i,],
                                             prob=p.hat[i,],
                                             size=size[i,],
                                             theta=exp(logtheta),
                                             log=TRUE))
    if (se) {
      o <- optim(par=1, fn=f, i=i, method="L-BFGS-B", lower=minld, upper=maxld,
                 control=list(fnscale=-1), hessian=TRUE)
      theta.hat[i] <- exp(o$par)
      var.est <- -1 * o$hessian^-1
      se.vec[i] <- if (var.est >= 0) sqrt(var.est) else NA
    } else {
      o <- optimize(f, interval=c(minld,maxld), i=i, maximum=TRUE)
      theta.hat[i] <- exp(o$maximum)
    }
  }
  if (se) {
    return(cbind(theta=theta.hat, se=se.vec))
  }
  theta.hat
}
