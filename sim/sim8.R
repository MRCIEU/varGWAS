library("dplyr")
library("broom")
source("funs.R")
set.seed(134)

n_sim <- 200
n_obs <- 1000

results <- data.frame()
for (i in 1:n_sim){
    # simulate covariates
    data <- data.frame(
        S = paste0("S", seq(1, n_obs)),
        X = get_simulated_genotypes(af, n_obs),
        stringsAsFactors=F
    )

    # simulate outcome
    data$Y <- rnorm(n_obs, sd=1 + data$X * 1)
    
    # run models
    res <- run_models(data)

    results <- rbind(results, res)
}