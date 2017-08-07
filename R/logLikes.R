#' Log likelihood for Negative Binomial
#'
#' This is a simple function to be passed to \code{apeglm}
#' as a log likelihood for the negative binomial distribution.
#'
#' @param y the counts
#' @param x a design matrix
#' @param beta the coefficient vector (natural log scale)
#' @param param the overdispersion (alpha in DESeq2 notation)
#' @param offset an offset matrix (natural log scale)
#'
#' @return the log likelihood for each sample in \code{y}
#'
#' @importFrom stats dnbinom
#' @export
#' 
#' @examples 
#' 
#' # this function is used by 'apeglm' to specify
#' # a negative binomial log likelihood.
#' # so it's only passed as an argument, not for use on its own.
#' # we can show its output nevertheless:
#'
#' y <- rnbinom(10, mu=100, size=1/.1)
#' x <- cbind(rep(1,10),rep(0:1,each=5))
#' beta <- c(log(100),0)
#' param <- .1
#' offset <- rep(0, 10)
#' logLikNB(y, x, beta, param, offset)
#' 
logLikNB <- function(y, x, beta, param, offset) {
  xbeta <- x %*% beta + offset
  mean.hat <- exp(xbeta)
  dnbinom(y, mu=mean.hat, size=1/param, log=TRUE)
}
