# power simulation for UKBB data

library("data.table")
library("broom")
library("genpwr")
library('optparse')
library("jlst")
source("funs.R")
set.seed(12345)

n_obs <- 500000
n_sim <- 200
af <- 0.4

# main effect size of X on Y detectable with 80% power
delta <- 0.05

# interaction effect relative to main effect
theta <- delta * 1

# simulate GxE interaction effects and estimate power
results <- data.frame()  
for (i in 1:n_sim) {
    # simulate variables
    X <- get_simulated_genotypes(af, n_obs); X <- scale(X)
    U <- rnorm(n_obs)
    Y <- X * delta + U * delta + X * U * theta + rnorm(n_obs, sd=sqrt(1-(delta^2+delta^2+theta^2)))

    # run LAD-BF
    result <- vargwas_model(data.frame(X,Y), "X", "Y")
    result$varY <- var(Y)

    # store result
    results <- rbind(results, result)
}

# estimate power
binom.test(sum(results$phi_p < 0.05), n_sim) %>% tidy