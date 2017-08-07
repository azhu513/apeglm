#set.seed(100)
n.per.group <- 10 # number samples per group
n <- n.per.group * 4
m <- 1000 # number of rows

x <- factor(rep(letters[1:2],each=n.per.group*2))
z <- factor(rep(rep(letters[1:2],each=n.per.group),2))
design <- model.matrix(~ x + z + x:z) 

beta.sd <- 2
beta.x <- rnorm(m, 0, beta.sd)
beta.z <- rnorm(m, 0, beta.sd)
interaction <- rnorm(m, 0, beta.sd/2) # smaller distn
intercept <- 1
beta.mat <- cbind(intercept, beta.x, beta.z, interaction)

prob <- 1/(1 + exp(-t(design %*% t(beta.mat))))
size <- 10
Y <- matrix(rbinom(m*n, size=10, prob=prob),ncol=n)

# check if makes sense:
#one.row.y <- cbind(Y[1,],size-Y[1,])
#coef(glm(one.row.y ~ x + z + x:z, family=binomial))
#beta.mat[1,]
