library(devtools)
load_all()
## The first example - Binomial Data 
source("simulate_binom.R")
## define the log likelihood function
# In this function, y is read counts, x is the design matrix
# beta is the regression coefficients vector
# param is the size to be passed to dbinom()
my.log.lik <- function(y, x, beta, param, offset) {
  xbeta <- x %*% beta
  mean.hat <- exp(xbeta)/(1+exp(xbeta))
  dbinom(y, size = param, prob = mean.hat, log=TRUE)
}
param <- rep(10,nrow(Y))
#C <- matrix(c(0,0,0,1),nrow=4,ncol=1)
#colnames(C) <- "interaction"
alpha <- .01
source("simulate_binom.R")
system.time({
  fit <- apeglm(Y=Y, x=design, log.lik=my.log.lik,
                param=param, coef=4,
                interval.level=1-alpha,
                log.link=FALSE)
})
mean(beta.mat[,4] > fit$interval[,1] & beta.mat[,4] < fit$interval[,2], na.rm=TRUE)
lap.int <- matrix(fit$map[,4] + qnorm(rep(c(alpha/2,1-alpha/2),each=nrow(fit$map))) * fit$se[,4],ncol=2)
mean(beta.mat[,4] > lap.int[,1] & beta.mat[,4] < lap.int[,2])
sum(is.na(fit$interval[,1]))

library(rafalib)
bigpar(2,2)
for (i in seq_len(ncol(beta.mat))) {
  plot(beta.mat[,i], fit$map[,i])
}


#########################################################

## The second example - Negative binomial Data

source("simulate_NB.R")
load_all()

## define the log likelihood function
# In this function, y is read counts, x is the design matrix
# beta is the regression coefficients vector
# param is the dispersion parameter to be passed to dnbinom()
my.log.lik <- function(y, x, beta, param, offset) {
  if (is.null(offset)){
    xbeta <- x %*% beta
  } else {
    xbeta <- x %*% beta + offset
  }
  mean.hat <- exp(xbeta)
  dnbinom(y, mu=mean.hat, size=1/param, log=TRUE)
}

n.per.group <- 3
z <- sim(1,1000,n.per.group*2,.01)
Y <- z[[1]]
beta.mat <- z[[2]]
param <- rep(z[[3]], nrow(Y))
condition <- factor(rep(letters[1:2],each=n.per.group))
design <- model.matrix(~ condition)

coef <- 2
mean(beta.mat[,coef] > fit$interval[,1] & beta.mat[,coef] < fit$interval[,2], na.rm=TRUE)
lap.int <- matrix(fit$map[,coef] + qnorm(rep(c(alpha/2,1-alpha/2),each=nrow(fit$map))) * fit$se[,coef],ncol=2)
mean(beta.mat[,coef] > lap.int[,1] & beta.mat[,coef] < lap.int[,2])

#

plot(sort(fit$fsr))
svals <- svalue(fit$fsr)
prop.table(table(sign(fit$map[svals < .01,2]) != sign(beta.mat[svals < .01,2])))

#

plot(sort(-log10(fit$fsr)))

#

plot(fit$map[,coef]/fit$se[,coef], -log10(fit$fsr))
plot(beta.mat[,coef], fit$map[,coef])

#

hist(fit$map[,2],breaks=30,col="grey")
library(rafalib)
bigpar(2,2)
for (i in seq_len(ncol(beta.mat))) {
  plot(beta.mat[,i], fit$map[,i])
}

#

bigpar()
z <- fit$map[,coef] / fit$se[,coef]
p <- pnorm(-abs(z))
plot(fit$map[,coef]/fit$se[,coef], fit$fsr)
plot(p, fit$fsr, col=ifelse(z>0,"blue","red"))
abline(0,1);abline(v=0.5,h=0.5)
plot(-log10(p), -log10(fit$fsr), col=ifelse(z>0,"blue","red"))
abline(0,1)

#

s <- seq_len(nrow(Y))
plot((fit$interval[,1] - fit$map[,2]) / fit$se[,2], s, xlim=c(-5,5))
points((fit$interval[,2] - fit$map[,2]) / fit$se[,2], s)

#

s <- seq_len(nrow(Y))
o <- order(fit$map[,2])
plot(fit$map[o,2] + fit$se[o,2], s, xlim=c(-1,1))
points(fit$map[o,2] - fit$se[o,2], s, col="red")

###

# testing narrow distribution of betas
# how this affects FSR and MAE

pc <- list(
  noshrink = 1,
  prior.mean = 0,
  prior.scale = 1,
  prior.df = 1,
  prior.noshrink.mean = 0,
  prior.noshrink.scale = 15
)
pc0 <- pc

library(pbapply)
library(ggplot2)
alpha <- .01

m <- 500 # number of rows
dispersion <- .0001
s <- c(.01,.02,.03,.05,.1,.2,.3,.5,1)
mult <- c(1,1.5,2,2.5)
ns <- c(50,100) # n=3 not much error, n=100 big error
grid <- expand.grid(ns=ns,s=s,mult=mult)
grid <- data.frame(lapply(grid, rep, each=3))
res <- pbsapply(seq_len(nrow(grid)), function(i) {
  n.per.group <- grid$ns[i]; n <- n.per.group * 2
  condition <- factor(rep(letters[1:2],each=n.per.group))
  design <- model.matrix(~ condition)
  dat <- sim(grid$s[i],m,n,dispersion); Y <- dat[[1]]
  param <- rep(dat[[3]], nrow(Y))
  pc$prior.scale <- grid$mult[i] * grid$s[i]
  fit <- apeglm(Y=Y,x=design,log.lik=my.log.lik,param=param,prior.control=pc,coef=2)
  fit2 <- apeglm(Y=Y,x=design,log.lik=my.log.lik,param=param,prior.control=pc0,coef=2)
  cols <- ifelse(dat[[2]][,1] < 5, "red", "black")
  plot(dat[[2]][,2], fit$map[,2], col=cols);abline(0,1,col="red")
  plot(dat[[2]][,2], fit2$map[,2], col=cols);abline(0,1,col="red")
  eif <- median(abs(dat[[2]][,2] - fit$map[,2]))/median(abs(dat[[2]][,2] - fit2$map[,2]))
  fsr <- if (sum(fit$svalue < alpha) == 0) { 0 } else {
           mean(sign(fit$map[fit$svalue < alpha,2]) !=
                sign(dat[[2]][fit$svalue < alpha,2]), na.rm=TRUE)
         }
  fsr2 <- if (sum(fit2$svalue < alpha) == 0) { 0 } else {
            mean(sign(fit2$map[fit2$svalue < alpha,2]) !=
                 sign(dat[[2]][fit2$svalue < alpha,2]), na.rm=TRUE)
          }
  c(eif, fsr, fsr2)
}, cl=2)
res <- t(res); colnames(res) <- c("eif","fsr","fsr2")
gres <- data.frame(grid, res)
gres$mult <- factor(gres$mult)
ggplot(gres, aes(s, fsr, col=mult, group=mult)) +
  stat_summary(fun.y=mean,fun.ymax=max,fun.ymin=min) +
  geom_smooth(se=FALSE,size=.5) +
  scale_x_log10(breaks=s) + facet_wrap(~ns) +
  geom_hline(yintercept=alpha, col="orange")
ggplot(gres, aes(s, eif, col=mult, group=mult)) +
  stat_summary(fun.y=mean,fun.ymax=max,fun.ymin=min) +
  geom_smooth(se=FALSE,size=.5) +
  scale_x_log10(breaks=s) + facet_wrap(~ns) +
  geom_hline(yintercept=1, col="orange")
ggplot(gres, aes(s, fsr2)) + geom_point(size=.75) +
  geom_smooth(se=FALSE,size=.5) +
  scale_x_log10(breaks=s) + facet_wrap(~ns) +
  geom_hline(yintercept=alpha, col="orange")
###

# JS estimator for Normal scale

Atrue <- .2
n <- 1000
res <- t(sapply(1:100, function(i) {
  set.seed(i)
  theta <- rnorm(n, 0, sqrt(Atrue))
  D <- rexp(n,rate=.1)
  X <- rnorm(n, theta, sqrt(D))
  (rough <- var(X) - mean(D))
  S <- X^2
  I <- function(A) 1/(2*(A + D)^2)
  Ahat <- function(A) sum((S - D) * I(A)) / sum(I(A))
  objective <- function(A) Ahat(A) - A
  minA <- .01
  if (objective(minA) < 0) {
    zero <- minA
  } else {
    zero <- uniroot(objective, interval=c(.01, 20))$root
  }
  #s <- seq(from=.01, to=2, by=.01)
  #plot(s,sapply(s, objective))
  c(rough, zero)
}))
head(res)
boxplot(res)
abline(h=Atrue, col="blue")

