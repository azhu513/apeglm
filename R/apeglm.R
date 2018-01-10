#' Approximate posterior estimation for GLM coefficients
#'
#' apeglm provides Bayesian shrinkage estimators for effect sizes
#' in GLM models, using approximation of the posterior for individual coefficients.
#'
#' \code{prior.control} is a list of parameters that will be passed to determine
#' the prior distribution. Users are allowed to have a Normal prior on the 
#' intercept, and a t prior on the non-intercept coefficients (similar
#' to \code{bayesglm} in the \code{arm} package. The following are defaults:
#'
#' \itemize{
#'   \item \code{no.shrink = 1}: index of the coefficient(s) not to shrink 
#'   \item \code{prior.mean = 0}: mean of t prior
#'   \item \code{prior.scale = 1}: scale of t prior
#'   \item \code{prior.df = 1}: df of t prior
#'   \item \code{prior.no.shrink.mean = 0}: mean of Normal
#'   \item \code{prior.no.shrink.scale = 15}: scale of Normal
#' }
#'
#' So without specifying \code{prior.control}, the following is set inside \code{apeglm}:
#'
#' \code{prior.control <- list(no.shrink=1,prior.mean=0,prior.scale=1,
#'       prior.df=1,prior.no.shrink.mean=0,prior.no.shrink.scale=15)}
#'
#' Note that the prior should be defined on the natural log scale for a log link GLM.
#' 
#' @param Y the observations, which can be a matrix or SummarizedExperiment,
#' with columns for samples and rows for "features" (e.g. genes in a genomic context).
#' If Y is a SummarizedExperiment, \code{apeglm} will return, in addition
#' to other list items, a GRanges or GRangesList \code{ranges} with the
#' estimated coefficients as metadata columns.
#' @param x design matrix, with intercept in the first column
#' @param log.lik the log of likelihood function, specified by the user.
#' For Negative Binomial distribution, user can use \code{logLikNB} provided within the package.
#' @param param the other parameter(s) to be used in the likelihood function,
#' e.g. the dispersion parameter for a negative binomial distribution.
#' this can be a vector or a matrix (with columns as parameters)
#' @param coef (optional) the index of the coefficient for which
#' to generate posterior estimates (FSR and intervals)
#' @param mle (optional) a 2 column matrix giving the MLE and its standard error
#' of \code{coef}. this will be used to adapt the scale of the prior (empirical Bayes).
#' This overrides the \code{prior.scale} specified by \code{prior.control}
#' and sets \code{no.shrink} to all coefficients other than \code{coef}.
#' Note that these MLE's and SE's should be on the natural log scale for a log link GLM.
#' @param no.shrink logical, if TRUE, apeglm won't perform shrinkage (default is FALSE)
#' @param interval.type (optional) can be "laplace", "HPD", or "credible", which specifies 
#' the type of Bayesian interval that the user wants to output; "laplace" represents the 
#' Laplace approximation of the posterior mode
#' @param interval.level (optional) default is 0.95
#' @param threshold (optional) a threshold for integrating posterior probabilities,
#' see details under 'Value'.
#' Note that this should be on the natural log scale for a log link GLM.
#' @param contrasts (optional) contrast matrix, same number of rows as \code{x}
#' @param weights (optional) weights matrix, same shape as \code{Y}
#' @param offset (optional) offsets matrix, same shape as \code{Y}.
#' Note that this should be on the natural log scale for a log link GLM.
#' @param flip.sign whether to flip the sign of threshold value
#' when MAP is negative, default is TRUE (threshold must then be positive)
#' @param prior.control see Details
#' @param multiplier a positive number, when the prior is adapted to the \code{mle}
#' matrix provided, this parameter connects the scale of the estimated distribution
#' of betas to the scale of the prior. the default value was chosen based on
#' FSR and error analysis of simulated data
#' @param ngrid the number of grid points for grid integration of intervals
#' @param nsd the number of SDs of the Laplace posterior approximation to set the
#' left and right edges of the grid around the MAP
#' @param ngrid.nuis the number of grid points for nuisance parameters
#' @param nsd.nuis the number of Laplace standard errors to set the
#' left and right edges of the grid around the MAP of the nuisance parameters
#' @param log.link whether the GLM has a log link (default = TRUE)
#' @param param.sd (optional) potential uncertainty measure on the parameter \code{param}.
#' this should only be a vector, used when \code{param} specifies a single parameter
#' @param method options for how apeglm will find the posterior mode and SD.
#' The default is "general" which allows the user to specify a likelihood
#' in a general way. Alternatives for faster performance with the Negative Binomial
#' likelihood are "negbinR" and "negbinC", which should provide increasing
#' benefits respectively in terms of speed. These second two options will
#' ignore any function provided to the \code{log.lik} argument, and \code{param}
#' should specify the dispersion parameter (such that Var = mu + param mu^2).
#' @param optim.method the method passed to \code{optim}
#' @param bounds the bounds for the numeric optimization 
#'  
#' @return a list of matrices containing the following components:
#' \itemize{
#'   \item \code{map}: matrix of MAP estimates, columns for coefficients and rows for features
#'   \item \code{sd}: matrix of posterior SD, same shape as \code{map}
#'   \item \code{prior.control}: list with details on the prior
#'   \item \code{fsr}: vector of the false sign rate for \code{coef}
#'   \item \code{interval}: matrix of either HPD or credible interval for \code{coef}
#'   \item \code{thresh}: vector of the posterior probability that the estimated parameter 
#' is smaller than the threshold value specified in \code{threshold}
#' when MAP is positive (or greater than
#' -1 * threshold value when MAP is negative and flip.sign is TRUE)
#'   \item \code{diag}: matrix of diagnostics
#'   \item \code{contrast.map}: vector of MAP estimates corresponding to the \code{contrast}
#' when \code{contrast} is given 
#'   \item \code{contrast.sd}: vector of posterior SD corresponding to the \code{contrast}
#' when \code{contrast} is given
#'   \item \code{ranges}: GRanges or GRangesList with the estimated coefficients,
#' if \code{Y} was a SummarizedExperiment.
#' }
#'
#' Note that all parameters associated with coefficients,
#' e.g. \code{map}, \code{sd}, etc., are returned on the natural log scale for a log link GLM.
#' 
#' @importFrom SummarizedExperiment assay rowRanges
#' @importFrom GenomicRanges mcols<-
#' @importFrom stats dnorm pnorm qnorm sd dt optim uniroot
#' @importFrom utils head tail
#' @importFrom methods is
#' 
#' @export
#' 
#' @examples 
#' 
#' # Simulate RNA-Seq read counts data
#' 
#' # 5 samples for each of the two groups
#' # a total of 100 genes
#' n.per.group <- 5 
#' n <- n.per.group * 2
#' m <- 100
#' 
#' # The design matrix includes one column of intercept
#' # and one column indicating samples that belong to the second group
#' condition <- factor(rep(letters[1:2], each = n.per.group))
#' x <- model.matrix(~condition) 
#' 
#' # Specify the standard deviation of beta (LFC between groups)
#' beta.sd <- 2
#' beta.cond <- rnorm(m, 0, beta.sd)
#' beta.intercept <- runif(m, 2, 6)
#' beta.mat <- cbind(beta.intercept, beta.cond)
#' 
#' # Generate the read counts
#' mu <- exp(t(x %*% t(beta.mat)))
#' Y <- matrix(rnbinom(m*n, mu=mu, size=1/.1), ncol = n)
#' 
#' # Here we will use the negative binomial log likelihood
#' # which is an exported function. See 'logLikNB' for details.
#' # For the NB:
#' # 'param' is the dispersion estimate (1/size)
#' # 'offset' can be used to adjust for size factors (log of size factors)
#' param <- matrix(0.1, nrow = m, ncol = 1)
#' offset <- matrix(0, nrow = m, ncol = n)
#' 
#' # Shrinkage estimator of betas:
#' # (for adaptive shrinkage, 'apeglm' requires 'mle' coefficients
#' # estimated with another software, or by first running 'apeglm'
#' # setting 'no.shrink=TRUE'.)
#' res <- apeglm(Y = Y, x = x,
#'               log.lik = logLikNB,
#'               param = param,
#'               offset = offset,
#'               coef = 2)
#'
#' head(res$map)
#' plot(beta.mat[,2], res$map[,2])
#' abline(0,1)
#' 
apeglm <- function(Y, x, log.lik, 
                   param=NULL,
                   coef=NULL,
                   mle=NULL,
                   no.shrink=FALSE,
                   interval.type=c("laplace", "HPD", "credible"),
                   interval.level=0.95,
                   threshold=NULL, contrasts,
                   weights=NULL, offset=NULL,
                   flip.sign=TRUE,
                   prior.control,
                   multiplier=1,
                   ngrid=50, nsd=5,
                   ngrid.nuis=5, nsd.nuis=2,
                   log.link=TRUE,
                   param.sd=NULL,
                   method=c("general","negbinR","negbinC"),
                   optim.method="BFGS",
                   bounds=c(-Inf,Inf)) {

  if (missing(prior.control)) {
    prior.control <- list(
      no.shrink = 1,
      prior.mean = 0,
      prior.scale = 1,
      prior.df = 1,
      prior.no.shrink.mean = 0,
      prior.no.shrink.scale = 15
    )
  }

  if (no.shrink) {
    prior.control$no.shrink <- seq_len(ncol(x))
  }

  stopifnot(multiplier > 0)
  
  if (!is.null(mle)) {
    stopifnot(!is.null(coef))
    prior.control$no.shrink <- setdiff(seq_len(ncol(x)), coef)
    prior.control$prior.var <- priorVar(mle)
    prior.scale <- multiplier * sqrt(prior.control$prior.var)
    prior.scale <- min(prior.scale, 1)
    prior.control$prior.scale <- prior.scale
  }

  stopifnot(ncol(Y) == nrow(x))
  interval.type <- match.arg(interval.type)
  method <- match.arg(method)
  
  if (!is.matrix(param)) param <- as.matrix(param, ncol=1)
  # don't have code yet for use of threshold with param.sd
  stopifnot(is.null(param.sd) | is.null(threshold))
  if (flip.sign==TRUE & !is.null(threshold)) stopifnot(threshold > 0)
  
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2]]
  
  nvars <- ncol(x)

  hasRanges <- FALSE
  if (is(Y, "SummarizedExperiment")) {
    if (is(Y, "RangedSummarizedExperiment")) {
      ranges <- rowRanges(Y)
      hasRanges <- TRUE
    }
    Y <- assay(Y)
  }
  
  Y <- as.matrix(Y)
  rownames <- dimnames(Y)[[1]]
  nrows <- nrow(Y)
  
  if (!is.null(weights)) {
    stopifnot(ncol(weights)==ncol(Y))
    stopifnot(nrow(weights)==nrow(Y))
  }

  if (!is.null(offset)) {
    stopifnot(ncol(offset)==ncol(Y))
    stopifnot(nrow(offset)==nrow(Y))
  }

  if (method == "general") {
    offset.in.log.lik <- any(grepl("offset",as.character(body(log.lik))))
    if (offset.in.log.lik) {
      if (is.null(offset)) {
        stop("log.lik uses 'offset', so 'offset' should be non-NULL")
      }
    }
  }
  
  intercept.idx <- rowSums(x == 0) == nvars - 1
  if (nrows >= 2) {
    intercept <- rowMeans(Y[,intercept.idx,drop=FALSE])
  } else {
    intercept <- mean(Y[,intercept.idx,drop=FALSE])
  }
  
  result <- list()
  result$map <- matrix(nrow=nrows, ncol=nvars,
                         dimnames=list(rownames, xnames))
  result$sd <- matrix(nrow=nrows, ncol=nvars,
                      dimnames=list(rownames, xnames))
  result$prior.control <- prior.control
  
  if (!is.null(coef)) {
    if (!(is.numeric(coef) & coef==round(coef) & length(coef) == 1)){
      stop("coef must be numeric vector of length 1, and an integer")
    }
    if (coef < 2 | coef > ncol(x)){
      stop("'coef' must be between 2 and the number of columns of 'x'")
    }
    result$fsr <- matrix(nrow=nrows, ncol=1, dimnames=list(rownames, xnames[coef]))
    result$svalue <- matrix(nrow=nrows, ncol=1, dimnames=list(rownames, xnames[coef]))
    interval.nms <- list(rownames, c(paste0((1-interval.level)/2*100,"%"),
                                     paste0((1-(1-interval.level)/2)*100,"%")))
    result$interval <- matrix(nrow=nrows, ncol=2, dimnames=interval.nms)
    if (!is.null(threshold)) {
      result$thresh <- matrix(nrow=nrows, ncol=1,
                              dimnames=list(rownames, xnames[coef]))
    }
    if (interval.type != "laplace") {
      result$diag <- matrix(NA, nrow=nrows, ncol=4, 
                            dimnames=list(rownames,
                                          c("conv","count","out.left","out.right")))
    } else {
      result$diag <- matrix(NA, nrow=nrows, ncol=2,
                            dimnames=list(rownames, c("conv","count")))
    }
  } else {
    result$diag <- matrix(NA, nrow=nrows, ncol=2,
                          dimnames=list(rownames, c("conv","count")))
  }
  if (!missing(contrasts)) {
    ncontr <- ncol(contrasts)
    contrast.nms <- list(rownames, colnames(contrasts))
    result$contrast.map <- matrix(nrow=nrows, ncol=ncontr,
                                  dimnames=contrast.nms)
    result$contrast.sd <- matrix(nrow=nrows, ncol=ncontr,
                                 dimnames=contrast.nms)
  }
    
  for (i in seq_len(nrows)) {
    weights.row <- if (is.null(weights)) NULL else weights[i,]
    offset.row <- if (is.null(offset)) NULL else offset[i,]
    param.i <- if (is.null(param)) NULL else param[i,,drop=TRUE] # drop the dimension
    param.sd.i <- if (is.null(param.sd)) NULL else param.sd[i]
    row.result <- apeglm.single(y = Y[i,], x=x, log.lik=log.lik, 
                                param=param.i,
                                coef=coef,
                                interval.type=interval.type,
                                interval.level=interval.level,
                                threshold=threshold,
                                contrasts=contrasts, 
                                weights=weights.row, offset=offset.row,
                                flip.sign = flip.sign,
                                prior.control=prior.control,
                                ngrid=ngrid, nsd=nsd,
                                ngrid.nuis=ngrid.nuis, nsd.nuis=nsd.nuis,
                                log.link=log.link,
                                param.sd=param.sd.i,
                                intercept=intercept[i],
                                method=method,
                                optim.method=optim.method,
                                bounds=bounds)
    result$map[i,] <- row.result$map
    result$sd[i,] <- row.result$sd
    if (!is.null(coef)) {
      result$fsr[i,] <- row.result$fsr
      result$interval[i,] <- row.result$ci
      if (!is.null(threshold)){
        if (flip.sign == TRUE & row.result$map[coef] < 0){ 
          result$thresh[i,] <- 1 - row.result$thresh
        } else {
          result$thresh[i,] <- row.result$thresh
        }
      }
    }
    result$diag[i,] <- row.result$diag
    if (!missing(contrasts)) { 
      result$contrast.map[i,] <- row.result$contrast.map
      result$contrast.sd[i,] <- row.result$contrast.sd
    }
  }

  if (!is.null(coef)) {
    result$svalue[,1] <- svalue(result$fsr)
  }

  if (hasRanges) {
    mcols(ranges) <- result$map
    result$ranges <- ranges
  }

  return(result)
}


# Log of prior density
# Allow for t or normal prior for intercept, and t for coefficients
# work on scale later
log.prior <- function(beta, prior.control) {
  p <- prior.control
  log.prior.no.shrink <- sum(dnorm(beta[p$no.shrink], mean = p$prior.no.shrink.mean, 
                                   sd = p$prior.no.shrink.scale, log=TRUE))
  # the scaling of the t density can be left out after taking the log
  log.prior.shrink <- sum(dt(beta[-p$no.shrink]/p$prior.scale,
                             df = p$prior.df, ncp = p$prior.mean, log=TRUE))
  log.prior.no.shrink + log.prior.shrink
}

log.post <- function(beta, log.lik, log.prior, y, x, param, weights, 
                     offset, prior.control) {
  
  if (is.null(weights)) {
    out <- sum(log.lik(y, x, beta, param, offset)) + log.prior(beta, prior.control)
  } else {
    out <- t(weights) %*% log.lik(y, x, beta, param, offset) + log.prior(beta, prior.control)
  }
  if (!is.finite(out)) return(-1e100)
  out
}


# For each row, we use this apeglm.single
# (unexported)
apeglm.single <- function(y, x, log.lik, 
                          param,
                          coef,
                          interval.type,
                          interval.level,
                          threshold, contrasts,
                          weights, offset,
                          flip.sign,
                          prior.control,
                          ngrid, nsd,
                          ngrid.nuis, nsd.nuis,
                          log.link=TRUE,
                          param.sd,
                          intercept,
                          method,
                          optim.method,
                          bounds) {

  if (log.link & all(y == 0)) {
    out <- buildNAOut(coef, interval.type, threshold, contrasts)
    return(out)
  }
    
  init <- numeric(ncol(x))
  if (log.link) {
    if (intercept == 0) {
      init[1] <- 0
    } else {
      init[1] <- log(intercept)
    }
  }

  # numerical optimization to find the MAP and posterior SD
  if (optim.method != "L-BFGS-B") {
    bounds <- c(-Inf, Inf)
  }

  if (method == "general") {
    o <- optim(par = init, fn = log.post, log.lik = log.lik, log.prior = log.prior, 
               y = y, x = x, param = param,
               weights = weights, offset = offset, 
               prior.control = prior.control, 
               control=list(fnscale=-1),
               lower=bounds[1], upper=bounds[2],
               hessian=TRUE, method=optim.method)
  } else if (method == "negbinR") {
    o <- optimNegBin(init=init, y=y, x=x, param=param,
                     weights=weights, offset=offset,
                     prior.control=prior.control,
                     bounds=bounds,
                     optim.method=optim.method)
  }
  
  map <- o$par
  sigma <- -solve(o$hessian)

  if (any(diag(sigma) <= 0)) {
    out <- buildNAOut(coef, interval.type, threshold, contrasts)
    out$map <- map
    return(out)
  }

  sd <- sqrt(diag(sigma))

  out <- list(map=map, sd=sd)
  # calculate statistics for a particular coefficient
  if (!is.null(coef)) {
    if (interval.type == "laplace") {
      stopifnot(is.null(param.sd)) # not implemented
      qn <- qnorm((1 - interval.level)/2,lower.tail=FALSE)
      out$ci <- c(map[coef] - qn * sd[coef], map[coef] + qn * sd[coef])
      out$diag <- c(o$convergence, o$counts[1])
      out$fsr <- pnorm(-abs(map[coef]),0,sd[coef])
      if (!is.null(threshold)) {
        if (flip.sign) {
          out$threshold <- pnorm(sign(map[coef])*threshold, map[coef], sd[coef])
        } else {
          out$threshold <- pnorm(threshold, map[coef], sd[coef])
        }
      }
    } else {
      # use a grid to evaluate posterior
      out <- gridResults(y=y, x=x, log.lik=log.lik, 
                         param=param, 
                         coef=coef,
                         interval.type=interval.type,
                         interval.level=interval.level,
                         threshold=threshold, contrasts=contrasts,
                         weights=weights, offset=offset,
                         flip.sign=flip.sign,
                         prior.control=prior.control,
                         ngrid=ngrid, nsd=nsd,
                         ngrid.nuis=ngrid.nuis, nsd.nuis=nsd.nuis,
                         log.link=log.link,
                         param.sd=param.sd,
                         intercept=intercept,
                         o=o, map=map, sigma=sigma, sd=sd, out=out)
    }
  } else {
    # just diagnostic values (if coef not specified)
    out$diag <- c(o$convergence, o$counts[1])
  }
  
  # calculate contrasts
  if (!missing(contrasts)) {
    contrasts <- data.matrix(contrasts)
    stopifnot(nrow(contrasts) == ncol(x))
    stopifnot(ncol(contrasts) >= 1)
    out$contrast.map <- map %*% contrasts
    out$contrast.sd <- t(sqrt(diag(t(contrasts) %*% sigma %*% contrasts)))
  }
  
  out
}

buildNAOut <- function(coef, interval.type, threshold, contrasts) {
  out <- list(map=NA, sd=NA)
  if (!is.null(coef)) {
    out$diag <- if (interval.type != "laplace") {
                  c(NA, NA, NA, NA)
                } else {
                  c(NA, NA)
                }
    out$fsr <- NA
    out$ci <- c(NA, NA)
    if (!is.null(threshold)) {
      out$threshold <- NA
    }
  } else {
    out$diag <- c(NA, NA)
  }
  if (!missing(contrasts)) {
    out$contrast.map <- NA
    out$contrast.sd <- NA
  }
  out
}

optimNegBin <- function(init, y, x, param, weights, offset, prior.control,
                        bounds, optim.method) {
  # TODO: weights and offset are being ignored
  size <- 1/param
  no.shrink <- prior.control$no.shrink
  shrink <- setdiff(seq_along(init), no.shrink)
  sigma <- prior.control$prior.no.shrink.scale
  S <- prior.control$prior.scale
  f <- function(beta, x, y, size, sigma, S, no.shrink, shrink, const) {
    xbeta <- x %*% beta
    prior <- sum(-beta[no.shrink]^2/(2*sigma^2) - log(1 + beta[shrink]^2/S^2))
    -sum(y * xbeta - (y + size) * log(size + exp(xbeta))) - prior + const
  }
  const <- -f(init, x, y, size, sigma, S, no.shrink, shrink, 0) - 1
  gr <- function(beta, x, y, size, sigma, S, no.shrink, shrink, const) {
    xbeta <- x %*% beta
    prior <- sum(-beta[no.shrink]/sigma^2 - 2*beta[shrink]/(S^2 + beta[shrink]^2))
    -t(x) %*% (y - (y + size) * exp(xbeta) / (size + exp(xbeta)))
  }
  optim(init, f, gr=gr, x=x, y=y, size=size,
        sigma=sigma, S=S, no.shrink=no.shrink,
        shrink=shrink, const=const,
        lower=bounds[1], upper=bounds[2],
        hessian=TRUE, method=optim.method)
}
