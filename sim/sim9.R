library("dplyr")
library("broom")
library("ggplot2")
source("funs.R")
set.seed(2445)

n_sim <- 200
n_obs <- 5000

results <- data.frame()
for (i in 1:n_sim){
    # simulate covariates
    data <- data.frame(
        S = paste0("S", seq(1, n_obs)),
        X = get_simulated_genotypes(.4, n_obs),
        stringsAsFactors=F
    )

    # simulate outcome
    data$Y <- rnorm(n_obs, sd=sqrt(2 + data$X * 2))
    
    # run models
    res <- run_osca(data, T)
    res$dp <- dummy_p(data$X, data$Y)
    res$xp <- xsq_p(data$X, data$Y)
    oe <- get_osca_effect(res$dp,.4,n_obs,-1)
    res$oe_beta <- oe[1]
    res$oe_se <- oe[2]

    results <- rbind(results, res)
}

# test for correlation between methods
cor.test(results$dp, results$xp)
cor.test(results$dp, results$P.osca_median)