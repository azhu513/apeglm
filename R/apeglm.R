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
#' likelihood are:
#' "nbinomR", "nbinomCR", and "nbinomC" / "nbinomC*"
#' These alternative methods should provide increasing speeds, respectively.
#' From testing on RNA-seq data, they are roughly 5x, 10x and 50x faster than "general".
#' Note that "nbinomC" uses C++ to find the MAP for the coefficients,
#' but does not calculate or return the posterior SD or other quantities.
#' "nbinomC*" is the same as "nbinomC", but includes a random start for finding the MAP.
#' "nbinomCR" uses C++ to calculate the MAP and then estimates
#' the posterior SD in R, with the exception that if the MAP from C++
#' did not converge or gives negative estimates of posterior variance,
#' then this row is refit using optimization in R.
#' These alternatives require the degrees of freedom for the prior distribution to be 1,
#' and will ignore any function provided to the \code{log.lik} argument.
#' \code{param} should specify the dispersion parameter of a Negative Binomial
#' (such that Var = mu + param mu^2).
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
#' @importFrom Rcpp sourceCpp
#' @importFrom SummarizedExperiment assay rowRanges
#' @importFrom GenomicRanges mcols<-
#' @importFrom stats dnorm pnorm qnorm sd dt optim optimHess uniroot
#' @importFrom utils head tail
#' @importFrom methods is
#'
#' @useDynLib apeglm
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
                   method=c("general","nbinomR","nbinomCR","nbinomC","nbinomC*"),
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

  interval.type <- match.arg(interval.type)
  method <- match.arg(method)
  
  stopifnot(ncol(Y) == nrow(x))
  stopifnot(multiplier > 0)  
  
  if (no.shrink) {
    prior.control$no.shrink <- seq_len(ncol(x))
  }
  
  if (!is.null(mle)) {
    stopifnot(!is.null(coef))
    prior.control$no.shrink <- setdiff(seq_len(ncol(x)), coef)
    prior.control$prior.var <- priorVar(mle)
    prior.scale <- multiplier * sqrt(prior.control$prior.var)
    prior.scale <- min(prior.scale, 1)
    prior.control$prior.scale <- prior.scale
  }

  if (method != "general") {
    stopifnot(prior.control$prior.df == 1)
    stopifnot(prior.control$prior.mean == 0)
    stopifnot(prior.control$prior.no.shrink.mean == 0)
  }
    
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
  basemean <- rowMeans(Y[,intercept.idx,drop=FALSE])
  
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
    result$fsr <- matrix(nrow=nrows, ncol=1,
                         dimnames=list(rownames, xnames[coef]))
    result$svalue <- matrix(nrow=nrows, ncol=1,
                            dimnames=list(rownames, xnames[coef]))
    interval.nms <- list(rownames, c(paste0((1-interval.level)/2*100,"%"),
                                     paste0((1-(1-interval.level)/2)*100,"%")))
    result$interval <- matrix(nrow=nrows, ncol=2, dimnames=interval.nms)
    if (!is.null(threshold)) {
      result$thresh <- matrix(nrow=nrows, ncol=1,
                              dimnames=list(rownames, xnames[coef]))
    }
  }

  # diagnostic columns differ if we use laplace or not for posterior estimation
  if (interval.type != "laplace") {
    diag.cols <- c("conv","count","value","out.left","out.right")
  } else {
    diag.cols <- c("conv","count","value")
  }
  result$diag <- matrix(NA, nrow=nrows, ncol=length(diag.cols), 
                        dimnames=list(rownames, diag.cols))

  if (!missing(contrasts)) {
    ncontr <- ncol(contrasts)
    contrast.nms <- list(rownames, colnames(contrasts))
    result$contrast.map <- matrix(nrow=nrows, ncol=ncontr,
                                  dimnames=contrast.nms)
    result$contrast.sd <- matrix(nrow=nrows, ncol=ncontr,
                                 dimnames=contrast.nms)
  }

  # TODO eventually, break this out in sub-function
  if (method %in% c("nbinomCR","nbinomC","nbinomC*")) {
    nonzero <- rowSums(Y) > 0
    # the C++ code uses transposed data (samples x genes)
    YNZ <- t(Y[nonzero,,drop=FALSE])
    if (is.null(weights)) {
      weights <- matrix(1, nrow=nrow(Y), ncol=ncol(Y))
    }
    if (is.null(offset)) {
      offset <- matrix(0, nrow=nrow(Y), ncol=ncol(Y))
    }
    weightsNZ <- t(weights[nonzero,,drop=FALSE])
    offsetNZ <- t(offset[nonzero,,drop=FALSE])
    size <- 1/param[nonzero]
    sigma <- prior.control$prior.no.shrink.scale
    S <- prior.control$prior.scale
    no.shrink <- prior.control$no.shrink
    shrink <- setdiff(seq_len(ncol(x)), no.shrink)
    # now, estimate the scale of the function
    init <- rep(0, ncol(x))
    cnst <- sapply(seq_len(sum(nonzero)), function(i) {
      nbinomFn(init, x=x, y=YNZ[,i], size=size[i], weights=weightsNZ[,i],
               offset=offsetNZ[,i], sigma=sigma, S=S, no.shrink=no.shrink,
               shrink=shrink, cnst=0)
    })
    cnst <- ifelse(cnst > 1, cnst, 1)
    # now optimize over all rows using L-BFGS run in C++
    # on a scaled version of the negative posterior.
    # we run it twice to check for stability and issues w/ local maxima
    if (method == "nbinomC*") {
      init <- rnorm(ncol(x),0,.5)
    } else {
      init <- rep(c(.1,-.1),length.out=ncol(x))
    }
    out <- nbinomGLM(x=x, Y=YNZ, size=size, weights=weightsNZ,
                     offset=offsetNZ, sigma2=sigma^2, S2=S^2,
                     no_shrink=no.shrink, shrink=shrink,
                     init=init, cnst=cnst)
    if (method == "nbinomCR") {
      init2 <- rep(c(-.1,.1),length.out=ncol(x))
      out2 <- nbinomGLM(x=x, Y=YNZ, size=size, weights=weightsNZ,
                        offset=offsetNZ, sigma2=sigma^2, S2=S^2,
                        no_shrink=no.shrink, shrink=shrink,
                        init=init2, cnst=cnst)
    }
    ## valueR <- sapply(seq_len(sum(nonzero)), function(i) {
    ##   nbinomFn(out$beta[,i], x=x, y=YNZ[,i], size=size[i], weights=weightsNZ[,i],
    ##            offset=offsetNZ[,i], sigma=sigma, S=S, no.shrink=no.shrink,
    ##            shrink=shrink, cnst=0)/cnst[i] + 10
    ## })
    ## nas <- rep(NA, nrow(result$diag)); result$diag <- cbind(result$diag, valueR=nas)
    ## result$diag[nonzero,"valueR"] <- valueR
    result$map[nonzero,] <- t(out$betas)
    result$diag[nonzero,"conv"] <- out$convergence
    result$diag[nonzero,"value"] <- out$value
    if (method == "nbinomCR") {
      # if the two fits above disagree by .01, say it did not converge
      delta <- apply(abs(out$betas - out2$betas), 2, max)
      result$diag[nonzero,"conv"][delta > .01] <- -1
    } else {
      # nbinomC or nbinomC* just return the result
      return(result)
    }
  }
  
  for (i in seq_len(nrows)) {
    weights.row <- if (is.null(weights)) NULL else weights[i,]
    offset.row <- if (is.null(offset)) NULL else offset[i,]
    param.i <- if (is.null(param)) NULL else param[i,,drop=TRUE] # drop the dimension
    param.sd.i <- if (is.null(param.sd)) NULL else param.sd[i]
    prefit.beta <- if (method == "nbinomCR") result$map[i,] else NULL
    prefit.conv <- if (method == "nbinomCR") result$diag[i,"conv"] else NULL

    row.result <- apeglm.single(y = Y[i,], x=x, log.lik=log.lik, 
      param=param.i, coef=coef, interval.type=interval.type, interval.level=interval.level,
      threshold=threshold, contrasts=contrasts, weights=weights.row, offset=offset.row,
      flip.sign=flip.sign, prior.control=prior.control,
      ngrid=ngrid, nsd=nsd, ngrid.nuis=ngrid.nuis, nsd.nuis=nsd.nuis,
      log.link=log.link, param.sd=param.sd.i, basemean=basemean[i],
      prefit.beta=prefit.beta, prefit.conv=prefit.conv,
      method=method, optim.method=optim.method, bounds=bounds)
    
    result$map[i,] <- row.result$map
    result$sd[i,] <- row.result$sd
    if (!is.null(coef)) {
      result$fsr[i,] <- row.result$fsr
      result$interval[i,] <- row.result$ci
      if (!is.null(threshold) & !is.na(row.result$map[coef])) {
        if (flip.sign == TRUE & row.result$map[coef] < 0){ 
          result$thresh[i,] <- 1 - row.result$threshold
        } else {
          result$thresh[i,] <- row.result$threshold
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
apeglm.single <- function(y, x, log.lik, param, coef, interval.type, interval.level,
                          threshold, contrasts, weights, offset, flip.sign, prior.control,
                          ngrid, nsd, ngrid.nuis, nsd.nuis, log.link=TRUE, param.sd,
                          basemean, prefit.beta, prefit.conv, method, optim.method, bounds) {

  if (log.link & all(y == 0)) {
    out <- buildNAOut(coef, interval.type, threshold, contrasts)
    return(out)
  }

  if (is.null(prefit.beta)) {
    init <- rep(c(1,-1),length.out=ncol(x))
    if (log.link) {
      if (basemean == 0) {
        init[1] <- 0
      } else {
        init[1] <- log(basemean)
      }
    }
  } else {
    init <- prefit.beta
  }

  # numerical optimization to find the MAP and posterior SD
  if (optim.method != "L-BFGS-B") {
    bounds <- c(-Inf, Inf)
  }

  if (method == "general") {
    o <- optim(par = init, fn = log.post,
               log.lik = log.lik, log.prior = log.prior, 
               y = y, x = x, param = param,
               weights = weights, offset = offset, 
               prior.control = prior.control, 
               control=list(fnscale=-1),
               lower=bounds[1], upper=bounds[2],
               hessian=TRUE, method=optim.method)
  } else if (method == "nbinomR") {
    # this optimizes all rows in R, with function and gradient
    # written specifically for negative binomial likelihood
    o <- optimNbinom(init=init, y=y, x=x, param=param,
                     weights=weights, offset=offset,
                     prior.control=prior.control,
                     bounds=bounds,
                     optim.method=optim.method)
  } else if (method == "nbinomCR") {
    # this uses previously estimated C++ MAP for negative binomial likelihood
    # and calculates the Hessian uses optimHess().
    # if the Hessian gives negative variance estimates for a row, it will
    # re-optimize the negative log posterior using method="nbinomR" for that row
    o <- optimNbinomHess(init=init, y=y, x=x, param=param,
                         weights=weights, offset=offset,
                         prior.control=prior.control,
                         bounds=bounds,
                         optim.method=optim.method,
                         prefit.conv=prefit.conv)
  }
  
  map <- o$par
  cov.mat <- -solve(o$hessian)

  if (any(diag(cov.mat) <= 0)) {
    out <- buildNAOut(coef, interval.type, threshold, contrasts)
    out$map <- map
    return(out)
  }

  sd <- sqrt(diag(cov.mat))

  out <- list(map=map, sd=sd)
  # calculate statistics for a particular coefficient
  if (!is.null(coef)) {
    # this is the default interval type - Laplace approximation of posterior
    if (interval.type == "laplace") {
      stopifnot(is.null(param.sd)) # not implemented
      qn <- qnorm((1 - interval.level)/2,lower.tail=FALSE)
      out$ci <- c(map[coef] - qn * sd[coef], map[coef] + qn * sd[coef])
      out$diag <- c(o$convergence, o$counts[1], o$value)
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
                         o=o, map=map, cov.mat=cov.mat, sd=sd, out=out)
    }
  } else {
    # just diagnostic values (if coef not specified)
    out$diag <- c(o$convergence, o$counts[1], o$value)
  }
  
  # calculate contrasts
  if (!missing(contrasts)) {
    contrasts <- data.matrix(contrasts)
    stopifnot(nrow(contrasts) == ncol(x))
    stopifnot(ncol(contrasts) >= 1)
    out$contrast.map <- map %*% contrasts
    out$contrast.sd <- t(sqrt(diag(t(contrasts) %*% cov.mat %*% contrasts)))
  }
  
  out
}

buildNAOut <- function(coef, interval.type, threshold, contrasts) {
  out <- list(map=NA, sd=NA)
  if (!is.null(coef)) {
    out$fsr <- NA
    out$ci <- c(NA, NA)
    if (!is.null(threshold)) {
      out$threshold <- NA
    }
  }
  if (interval.type != "laplace") {
    out$diag <- c(NA, NA, NA, NA, NA)
  } else {
    out$diag <- c(NA, NA, NA)
  }
  if (!missing(contrasts)) {
    out$contrast.map <- NA
    out$contrast.sd <- NA
  }
  out
}
