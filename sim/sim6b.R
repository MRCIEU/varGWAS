library("broom")
library("dplyr")
source("funs.R")
set.seed(23)

n_obs <- 1000
n_sim <- 200

results <- data.frame()
for (b in seq(0, 6, .5)){
    for (i in 1:n_sim){
        x <- get_simulated_genotypes(0.4, n_obs)
        u <- rnorm(n_obs)
        y <- x*u*b + rnorm(n_obs)

        # run models
        res <- run_models(data.frame(
            S = paste0("S", seq(1, n_obs)),
            X = x,
            U = u,
            Y = y,
            stringsAsFactors=F
        ))
        # estimate var(Y|G==1)
        b1 <- 1 * res$BETA_x.osca_median
        b2 <- 2 * res$BETA_x.osca_median
        se1 <- 1 * res$SE_x.osca_median
        se2 <- 2 * res$SE_x.osca_median
        lci1 <- b1 - (1.96 * se1)
        uci1 <- b1 + (1.96 * se1)
        lci2 <- b2 - (1.96 * se2)
        uci2 <- b2 + (1.96 * se2)
        # estimate var(Y|G==2)
        # store results
        results <- rbind(results, data.frame(
            b,
            t1=var(y[x==1]) - var(y[x==0]), t2=var(y[x==2]) - var(y[x==0]),
            b1, b2,
            lci1, uci1,
            lci2, uci2
        ))
    }
}

# write out results
write.table(file="osca-effects.txt", results)