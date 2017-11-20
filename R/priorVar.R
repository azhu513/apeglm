priorVar <- function(mle, min.var=.001^2, max.var=20^2) {
  # here estimate the variance of a Normal in a hierarchical model
  # from MLE estimates and their SE (so the non-equal variance hierarchical model)
  # here using notation and estimator from Efron and Morris:
  # "Data Analysis Using Stein's Estimator and its Generalizations" 1975
  stopifnot(ncol(mle) == 2) # expecting MLE estimates and SE as two columns of matrix
  keep <- !is.na(mle[,1])
  X <- mle[keep,1]
  D <- mle[keep,2]^2
  S <- X^2
  I <- function(A) 1/(2*(A + D)^2)
  Ahat <- function(A) sum((S - D) * I(A)) / sum(I(A))
  objective <- function(A) Ahat(A) - A
  if (objective(min.var) < 0) {
    zero <- min.var
  } else {
    # TODO this needs to have better error message on fail
    zero <- uniroot(objective, interval=c(min.var, max.var))$root
  }
  zero
}
