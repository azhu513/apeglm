#set.seed(1)
sim <- function(s, m, n, dispersion) {
  n.per.group <- n/2
  condition <- factor(rep(letters[1:2],each=n.per.group))
  design <- model.matrix(~ condition)  
  beta.sd <- s
  beta.b <- rnorm(m,0,beta.sd)
  #beta.c <- rnorm(m,0,beta.sd)
  intercept <- runif(m,1,10)
  beta.mat <- cbind(intercept, beta.b)
  #beta.mat <- cbind(intercept, beta.b, beta.c)
  mu <- exp(t(design %*% t(beta.mat)))
  Y <- matrix(rnbinom(m*n, size=1/dispersion, mu=mu),ncol=n)
  list(Y, beta.mat, dispersion)
}

# check if makes sense:
# library(MASS)
# coef(glm.nb(Y[1,] ~ condition))
# beta.mat[1,]
