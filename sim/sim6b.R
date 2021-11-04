library("broom")
library("dplyr")
source("funs.R")
set.seed(23)

n_obs <- 1000
n_sim <- 200

results <- data.frame()
for (b in seq(2,2)){
    for (i in 1:n_sim){
        # SNP
        x <- rbinom(n_obs, 2, .6)
        # modifier
        u <- rnorm(n_obs)
        # outcome
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
        e1 <- 1 * res$BETA_x.osca_median
        e2 <- 2 * res$BETA_x.osca_median
        se1 <- 1 * res$SE_x.osca_median
        se2 <- 2 * res$SE_x.osca_median
        lci1 <- e1 - (1.96 * se1)
        uci1 <- e1 + (1.96 * se1)
        lci2 <- e2 - (1.96 * se2)
        uci2 <- e2 + (1.96 * se2)
        # estimate var(Y|G==2)
        # store results
        results <- rbind(results, data.frame(
            v1=var(y[x==1]), # true variance of SNP=1 group
            v2=var(y[x==2]), # true variance of SNP=2 group
            e1, e2, se1, se2,
            lci1, uci1,
            lci2, uci2,
            b
        ))
    }
}

# check for coverage of CI
results %>% dplyr::group_by(b) %>%
    dplyr::summarize(tidy(binom.test(sum(v1 >= lci1 & v1 <= uci1), n()))) # count number of times variance of Y is within 95% CI

results %>% dplyr::group_by(b) %>%
    dplyr::summarize(tidy(binom.test(sum(v2 >= lci2 & v2 <= uci2), n()))) # count number of times variance of Y is within 95% CI

# write out results
write.table(file="osca-effects.txt", results)