library("broom")
set.seed(1234)
n_obs <- 1000

# simulate some data for testing
x <- rnorm(n_obs)
c1 <- rbinom(n_obs, 2, .5)
c2 <- rnorm(n_obs)
y <- 4 + 0.6*x + 2*c1 + 0.3*c2 + rnorm(n_obs)
id <- paste0("S", seq(1, n_obs))
d <- data.frame(id, x, c1, c2, y)

# fit
tidy(lm(y~x+c1+c2, data=d))

# write out
write.table(d, sep=",", quote=F, row.names=F, file="data.csv")

# outlier
n <- 100
x <- rbinom(n, 2, .5)
y <- x + rnorm(n)
write.csv(data.frame(x,y), row.names=F, quote=F, file="data-outlier.csv")