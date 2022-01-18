library("broom")
library("ggplot2")
library("dplyr")
source("funs.R")
set.seed(123)

n_sim <- 200
n_obs <- 1000
b <- 0.09 # main effect size to have 80% power

# determine power of main effect
#results <- data.frame()
#for (i in 1:n_sim){
#    x <- rnorm(n_obs)
#    y <- x*b + rnorm(n_obs)
#    res <- lm(y ~ x) %>% tidy %>% dplyr::filter(term=="x") %>% dplyr::select(p.value)
#    results <- rbind(results, res)
#}
#binom.test(sum(results$p.value<0.05),n_sim)

n_obs <- 100000
x <- rnorm(n_obs)
u <- rnorm(n_obs)
y <- x + u + x*u + rnorm(n_obs)
fit <- lm(y ~ x)
d <- resid(fit)^2
plot(x,d)