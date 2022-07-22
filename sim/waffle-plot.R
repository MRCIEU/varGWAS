library("ggplot2")
library("reshape")
set.seed(8)

n_obs <- 36
x <- rbinom(n_obs, 2, 0.25)
u <- rbinom(n_obs, 1, 0.5)
y <- x + x*u + rnorm(n_obs)
n <- 1:n_obs
df <- data.frame(x, u, y, n)

df <- expand.grid(x = 0:5, y = 0:5)
df$z <- y
# default is compatible with geom_tile()
ggplot(df, aes(x, y, fill = z)) + geom_raster()