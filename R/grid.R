gridResults <- function(y, x, log.lik, 
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
                        log.link,
                        param.sd,
                        o, map, cov.mat, sd, out) {

  cor.mat <- cov.mat / outer(sd, sd)
  
  # evaluate the unnormalized posterior for 'coef' over a grid 
  beta.min <- map[coef] - nsd*sd[coef]
  beta.max <- map[coef] + nsd*sd[coef]
  stopifnot(ngrid %% 2 == 0)
  betas <- seq(beta.min, beta.max, length.out=(ngrid + 1))

  # we integrate the unnormalized posterior over the other coefficients than 'coef'
  unpo <- integrateOther(cor.mat=cor.mat, coef=coef, betas=betas, map=map, sd=sd,
                         log.lik = log.lik, log.prior = log.prior, 
                         y = y, x = x, param = param, param.sd = param.sd,
                         prior.control = prior.control, 
                         weights = weights, offset = offset,
                         ngrid.nuis = ngrid.nuis, nsd.nuis = nsd.nuis)

  # check that unpo is monotonic on either side of map
  if (!all(diff(unpo[betas < map[coef]]) >= 0) 
      | !all(diff(unpo[betas > map[coef]]) <= 0)) {
    out$diag <- c(o$convergence, o$counts[1], o$value, NA, NA)
    out$fsr <- NA
    out$ci <- c(NA, NA)
    if (!is.null(threshold)) {
      out$threshold <- NA
    }
    return(out)
  }
    
  # integration within the grid
  unpo.grid.integral <- simpsons(unpo, ngrid, beta.min, beta.max)

  # estimate the area outside of the grid on left and right
  h <- (beta.max - beta.min)/ngrid # grid unit

  # check that the first grid points are not zero
  if (unpo[1] > 0 & unpo[2] > 0) {
    # fit a,b for tail = exp(a + bx)
    # just use the last two values and draw line going outward from grid
    left.b <- (log(unpo[2]) - log(unpo[1]))/h
    stopifnot(left.b > 0)
    left.a <- log(unpo[1]) - left.b * beta.min
    # integral of c * exp(b*x) = c/b exp(b*x)
    left.c <- exp(left.a)
    if (!is.finite(left.c) | !is.finite(exp(left.b * beta.min))){
      un.out.left <- 0
      un.out.zero.left <- 0
      if (!is.null(threshold)) {
        un.out.thresh.left <- 0
      }
    } else {
      # unnormalized area out of grid
      un.out.left <- left.c / left.b * exp(left.b * beta.min)
      # unnormalized area from 0 to -Inf and to Inf
      # (used to calculate out.area and FSR when the grid doesn't contain 0)
      un.out.zero.left <- left.c / left.b 
      if (!is.null(threshold)) {
        un.out.thresh.left <- left.c / left.b * exp(left.b * threshold)
      }
    }
  } else {
    un.out.left <- 0
    un.out.zero.left <- 0
    if (!is.null(threshold)) {
      un.out.thresh.left <- 0
    }
  }

  # check that last grid points are not zero
  if (unpo[ngrid] > 0 & unpo[ngrid+1] > 0) {
    right.b <- (log(unpo[ngrid+1]) - log(unpo[ngrid]))/h
    stopifnot(right.b < 0)
    right.a <- log(unpo[ngrid+1]) - right.b * beta.max
    right.c <- exp(right.a)
    if (!is.finite(right.c) | !is.finite(exp(right.b * beta.max))){
      un.out.right <- 0
      un.out.zero.right <- 0
    } else {
      un.out.right <- -right.c / right.b * exp(right.b * beta.max)
      un.out.zero.right <- -right.c / right.b
    }
  } else {
    un.out.right <- 0
    un.out.zero.right <- 0
  }
  
  # add outside area to grid integral
  stopifnot(!is.na(un.out.left) & !is.na(un.out.right))
  total.integral <- unpo.grid.integral + un.out.left + un.out.right
  
  # normalize posterior and areas out of grid
  norm.post <- unpo/total.integral
  out.left <- un.out.left/total.integral
  out.right <- un.out.right/total.integral

  # calculate the interval estimate:
  # first we check if the out-of-grid area is < (1 - interval.level) = alpha
  if ((out.left + out.right) < (1 - interval.level)) {
    # if so, we can use either of these two functions to find
    # interval boundaries within the grid (this is the desired routine)
    if (interval.type == "HPD"){
      # HPD interval:
      out$ci <- hpd.int(density = norm.post, grid = betas, h = h, 
                        level = interval.level)
    } else {
      # credible interval:
      out$ci <- cred.int(density = norm.post, grid = betas, h = h, 
                         level = interval.level, area.left = out.left)
    }
  } else {
    # if the out-of-grid area is > (1 - interval.level) = alpha
    # then we need to be more careful to find the left and/or right boundary
    #
    # we will estimate left and right tails to give us alpha/2 each.
    # if the area out-of-grid on left is less than alpha/2, this implies
    # we need to go into the grid to find the lower bound of credible interval.
    # otherwise we find the point out of the grid using an exponential
    # approximation of the tail.
    
    half.alpha <- (1 - interval.level)/2
    cumdensity <- cumsum(h*norm.post) + out.left
    if (out.left < half.alpha) {
      left.idx <- which(cumdensity > half.alpha)[1] - 1
      left.log.density <- log(norm.post[c(left.idx,(left.idx+1))])
      left.grid <- betas[c(left.idx,(left.idx+1))]
      intv.l <- solve.x(left.log.density, left.grid,
                        remn.area = half.alpha - cumdensity[left.idx])
    } else {
      intv.l <- solve.x(y = c(log(norm.post[2]), log(norm.post[1])),
                        x = c(betas[1], betas[2]),
                        remn.area = half.alpha)
    }
    if (out.right < half.alpha) {
      right.idx <- which(cumdensity > (1-half.alpha))[1] - 1
      right.log.density <- log(norm.post[c(right.idx,(right.idx+1))])
      right.grid <- betas[c(right.idx,(right.idx+1))]
      intv.u <- solve.x(right.log.density, right.grid,
                        remn.area = (1 - half.alpha) - cumdensity[right.idx])
    } else {
      intv.u <- solve.x(y = c(log(norm.post[ngrid]), log(norm.post[ngrid+1])),
                        x = c(betas[ngrid], betas[ngrid+1]),
                        remn.area = half.alpha)
    }
    out$ci <- c(intv.l, intv.u)
  }
  
  
  # diagnostic values including normalized tail area
  # to the left and right of the grid boundaries
  out$diag <- c(o$convergence, o$counts[1], o$value, out.left, out.right)
  
  # calculate the false sign rate
  
  # don't have code right now for FSR with param.sd
  if (!(is.null(param.sd))) {
    out$fsr <- NA
  } else {
    if (map[coef] > 0) {
      fsr.idx <- betas < 0
    } else {
      fsr.idx <- betas > 0
    }
    if (sum(fsr.idx) > 0 & is.null(param.sd)) {
      # see how much we over-integrate and remove this bit
      if (map[coef] > 0) {
        out.area <- out.left
        next.beta <- betas[max(which(fsr.idx)) + 1]
        over.spill.h <- next.beta - 0
        over.spill <- over.spill.h * tail(norm.post[fsr.idx],1)
      } else {
        out.area <- out.right
        next.beta <- betas[min(which(fsr.idx)) - 1]
        over.spill.h <- 0 - next.beta
        over.spill <- over.spill.h * head(norm.post[fsr.idx],1)
      }
      fsr <- sum(h*norm.post[fsr.idx]) - over.spill + out.area
    } else {
      # there are no beta grid points on the 'other side' of zero (relative to MAP)
      if (map[coef] > 0) {
        fsr <- un.out.zero.left/total.integral
      } else {
        fsr <- un.out.zero.right/total.integral
      }
    }

    out$fsr <- fsr
  }
  
  # provide posterior integral for a threshold
  # (other than 0 which is given by FSR)
  if (!is.null(threshold)) {
    if (flip.sign == TRUE & map[coef] < 0){
      threshold <- -threshold
    } 
    below.thresh.idx <- betas < threshold
    if (sum(below.thresh.idx) > 0) {
      next.beta <- betas[max(which(below.thresh.idx)) + 1]
      over.spill.h <- next.beta - threshold
      over.spill.thresh <- over.spill.h * tail(norm.post[below.thresh.idx],1)
      out$thresh <- sum(h*norm.post[below.thresh.idx]) - over.spill.thresh + out.left
    } else {
      out$thresh <- un.out.thresh.left/total.integral
    }
  }

  out
}

# Simpson's Rule 
# Intergration over a range
simpsons <- function(fx, ngrid, min, max) {
  h <- (max - min)/ngrid
  s <- fx[1] + fx[ngrid+1] + 2*sum(fx[seq(3,ngrid-1,by=2)]) + 4*sum(fx[seq(2,ngrid,by=2)])
  h/3 * s
}

cred.int <- function(density, grid, h, level, area.left) {
  
  prob <- (1-level)/2
  
  cumdensity <- cumsum(h*density) + area.left
  
  left.idx <- which(cumdensity > prob)[1] - 1
  right.idx <- which(cumdensity > (1-prob))[1] - 1
  
  left.log.density <- log(density[c(left.idx,(left.idx+1))])
  right.log.density <- log(density[c(right.idx,(right.idx+1))])
  
  left.grid <- grid[c(left.idx,(left.idx+1))]
  right.grid <- grid[c(right.idx,(right.idx+1))]
  
  intv.l <- solve.x(left.log.density, left.grid,
                    remn.area = prob - cumdensity[left.idx])
  intv.u <- solve.x(right.log.density, right.grid,
                    remn.area = (1 - prob) - cumdensity[right.idx])

  c(intv.l,intv.u)
}


## ref: "HDInterval"

hpd.int <- function(density, grid, h, level) {
  
  sorted.density <- sort(density, decreasing = TRUE)
  
  height.idx <- which(cumsum(sorted.density * h) >= level)[1]
  height <- sorted.density[height.idx]
  idcs <- which(density >= height)
  stopifnot(length(which(diff(idcs)>1))==0)
  
  half.remn <- (cumsum(sorted.density * h)[height.idx] - level)/2
  
  left.log.density <- log(density[idcs[1:2]])
  right.log.density <- log(density[idcs[c((length(idcs)-1),
                                          length(idcs))]])
  
  left.grid <- grid[idcs[1:2]]
  right.grid <- grid[idcs[c((length(idcs)-1),
                            length(idcs))]]
  
  intv.l <- solve.x(left.log.density, left.grid, remn.area = half.remn)
  intv.u <- solve.x(right.log.density, right.grid, remn.area = half.remn)
  
  c(intv.l,intv.u)
}


# function to find the point between x[1] and x[2]
# at which we achieve a given remaining area
# y is given in the log scale, but we are interested
# in the area between exp(y[1]) and exp(y[2])
solve.x <- function(y, x, remn.area) {
  b <- diff(y) / diff(x)
  c <- exp(y[1])
  # if we set the left grid point at x=0,
  # we want: int_0^x c exp(b*x) = remn.area
  # so c/b exp(b*x) - c/b = remn.area
  # then add this to x[1]
  x[1] + log( b/c * remn.area + 1)/b
}


integrateOther <- function(cor.mat, coef, betas, map, sd,
                           log.lik, log.prior, 
                           y, x, param, param.sd,
                           prior.control, 
                           weights, offset,
                           ngrid.nuis, nsd.nuis) {

  # skip the integration over other coefs if ngrid.nuis == 1
  if (ngrid.nuis == 1) {
    log.unpo <- sapply(betas, function(b) {
        beta.vec <- map
        beta.vec[coef] <- b
        log.post(beta.vec, log.lik = log.lik, log.prior = log.prior, 
                 y = y, x = x, param = param, 
                 prior.control = prior.control, 
                 weights = weights, offset = offset)
    })
    log.unpo <- log.unpo - max(log.unpo)
    unpo <- exp(log.unpo)
    return(unpo)
  }
  
  nnuis <- length(map) - 1
  nuis.grid.vec <- seq(-nsd.nuis, nsd.nuis, length.out=ngrid.nuis)
  nuis.grid <- as.matrix(expand.grid(lapply(seq_len(nnuis), function(i) nuis.grid.vec)))

  nuis.grid.big <- matrix(0, nrow=nrow(nuis.grid), ncol=length(map))
  nuis.grid.big[,-coef] <- nuis.grid

  ### old code using dmvnorm() to trim the number of grid evaluations ###
  ## if (nnuis > 1) {
  ##   nuis.probs <- dmvnorm(nuis.grid, rep(0,nnuis), cor.mat[-coef,-coef])
  ##   nuis.probs <- nuis.probs/sum(nuis.probs)
  ##   # only evaluate on the grid points larger than this cutoff
  ##   nuis.cutoff <- 1/nrow(nuis.grid)
  ##   stopifnot(any(nuis.probs > nuis.cutoff))
  ##   nuis.grid.big <- nuis.grid.big[nuis.probs > nuis.cutoff,]
  ## }

  nng <- nrow(nuis.grid.big)
  
  if (is.null(param.sd)) {
    log.unpo.mat <- matrix(nrow=nng, ncol=length(betas))
    log.unpo.mat <- fillLogUnpo(ng = nng,
                                LUM = log.unpo.mat,
                                betas = betas,
                                NGB = nuis.grid.big,
                                map = map, sd = sd, coef = coef,
                                log.lik = log.lik, log.prior = log.prior,
                                y = y, x = x, param = param,
                                prior.control = prior.control,
                                weights = weights, offset = offset)
    log.unpo.mat <- log.unpo.mat - max(log.unpo.mat)
    unpo.mat <- exp(log.unpo.mat)
    unpo <- colSums(unpo.mat)
  } else {
    # also include 'param' in the nuisance grid for integration
    # use a narrower grid than for the coefficients by 1/2 x
    npg <- ngrid.nuis
    param.grid <- seq(-nsd.nuis/2,nsd.nuis/2,length.out=npg)
    norm.po.mat <- matrix(nrow=npg, ncol=length(betas))
    for (j in seq_len(npg)) {
      param.star <- param + param.sd * param.grid[j]
      log.unpo.mat <- matrix(nrow=nng, ncol=length(betas))
      log.unpo.mat <- fillLogUnpo(ng = nng,
                                  LUM = log.unpo.mat,
                                  betas = betas,
                                  NGB = nuis.grid.big,
                                  map = map, sd = sd, coef = coef,
                                  log.lik = log.lik, log.prior = log.prior,
                                  y = y, x = x, param = param.star,
                                  prior.control = prior.control,
                                  weights = weights, offset = offset)
      log.unpo.mat <- log.unpo.mat - max(log.unpo.mat)
      unpo.mat <- exp(log.unpo.mat)
      unpo <- colSums(unpo.mat)
      unpo.grid.integral <- simpsons(unpo, length(betas)-1, min(betas), max(betas))
      # rough version of norm.post (ignores out of grid area)
      norm.po.mat[j,] <- unpo / unpo.grid.integral
    }
    # weight the contribution of the different evaluations
    # of 'param' by a normal density
    param.wts <- dnorm(param.grid)
    unpo <- colSums(norm.po.mat * param.wts)
  }
  
  unpo
}

# acronyms here
# LUM = log un-normalized posterior matrix
# NGB = nuisance grid, big
fillLogUnpo <- function(ng, LUM, betas, NGB, 
                        map, sd, coef, log.lik, log.prior, 
                        y, x, param , prior.control, 
                        weights, offset){
  for (i in seq_len(ng)) {
    LUM[i,] <- sapply(betas, function(b) {
      beta.vec <- map + NGB[i, ] * sd
      beta.vec[coef] <- b
      log.post(beta.vec, log.lik = log.lik, log.prior = log.prior, 
               y = y, x = x, param = param, 
               prior.control = prior.control, 
               weights = weights, offset = offset)})
  }
  LUM
}
