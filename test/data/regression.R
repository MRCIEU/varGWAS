library('quantreg')
set.seed(12345)

# create some data

n_obs <- 10000
mu <- 25
sigma <- 5
b <- 0.6

x <- rbinom(n_obs, 2, 0.5)
e <- rnorm(n_obs, mean = 0, sd = sigma)
y <- mu + b*x + e
df <- as.data.frame(cbind(x, y))

rqfit <- rq(y ~ x, data = df, tau = 0.5)
summary(rqfit)

write.csv(df, file="regression.csv", quote=F, row.names=F)