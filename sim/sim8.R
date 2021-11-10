library("dplyr")
library("broom")
source("funs.R")
set.seed(134)

n_sim <- 200
n_obs <- 10000

results <- data.frame()
for (i in 1:n_sim){
    # simulate covariates
    data <- data.frame(
        S = paste0("S", seq(1, n_obs)),
        X = get_simulated_genotypes(.4, n_obs),
        stringsAsFactors=F
    )

    # simulate outcome
    data$Y <- rnorm(n_obs, sd=sqrt(1 + data$X * 9))
    
    # run models
    res <- run_osca(data)
    res <- cbind(res, data.frame(
        v0=var(data$Y[data$X==0]),
        v1=var(data$Y[data$X==1]),
        v2=var(data$Y[data$X==2])
    ))
    bf <- dummy_model(data$X, data$Y)
    res$dummy1 <- bf[1]
    res$dummy2 <- bf[2]

    results <- rbind(results, res)
}