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

results <- data.frame()
for (i in 1:n_sim){
    for (phi in seq(0, 12, 0.5)){
        theta <- phi * b
        x <- get_simulated_genotypes(q, n_obs)
        u <- rnorm(n_obs)
        y <- x*b + x*u*theta + u*0 + rnorm(n_obs)
        results <- rbind(results, data.frame(
            phi, var=var(y) - 1, u="None"
        ))
        y <- x*b + x*u*theta + u*b + rnorm(n_obs)
        results <- rbind(results, data.frame(
            phi, var=var(y) - 1, u="Detectable with 80% power"
        ))
    }
}
