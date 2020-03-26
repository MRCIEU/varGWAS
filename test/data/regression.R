library('quantreg')
set.seed(12345)

# create some data
n_obs <- 50000
mu <- 25
sigma <- 5

# covariate-trait effects
b1 <- 0.6
b2 <- 2
b3 <- 0.05

# covariates
x <- rbinom(n_obs, 2, 0.5)
c1 <- rbinom(n_obs, 1, 0.5)
c2 <- runif(n_obs, min = 30, max = 70)

# create phenotype
e <- rnorm(n_obs, mean = 0, sd = sigma)
y <- mu + b1*x + b2*c1 + b3*c2 + e

# fit model
covar <- data.frame(c1=c1, c2=c2)
covar <- model.matrix(as.formula(~ .), data=covar)[,-1,drop=F]
fit <- suppressWarnings(rq(y~x+covar, tau=0.5))
summary(fit)
d <- abs(resid(fit))

# write out csv
df <- data.frame(x=x, c1=c1, c2=c2, y=y, d=d)
write.csv(df, file="regression.csv", quote=F, row.names=F)