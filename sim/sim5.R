library("dplyr")
library("broom")
library("tidyr")
library("ggpubr")
library("lmtest")
library("data.table")
source("funs.R")
set.seed(123)

# Requires OSCA and QCTOOL on PATH

n_sim <- 200
n_obs <- 1000

set_snp_confounder <- function(){
    r <- rep(NA, n_sim)
    p <- rep(NA, n_sim)
    for (i in 1:n_sim){
        C <- rnorm(n_obs)
        # confounder explains 5% variation of the SNP
        X <- sapply(C, function(c) {p <- 1 / (1 + exp(-c)); get_simulated_genotypes(p *.255, 1)})
        fit <- lm(X ~ C)
        r[i] <- summary(fit)$r.squared
        p[i] <- fit %>% tidy %>% dplyr::filter(term == "C") %>% dplyr::pull(p.value)
    }

    t.test(r)
    binom.test(sum(p<0.05), n_sim)
}

set_snp_outcome <- function(){
    px <- rep(NA, n_sim)
    rx <- rep(NA, n_sim)
    pc <- rep(NA, n_sim)
    rc <- rep(NA, n_sim)
    for (i in 1:n_sim){
        C <- rnorm(n_obs)
        X <- sapply(C, function(c) {p <- 1 / (1 + exp(-c)); get_simulated_genotypes(p *.255, 1)})
        # confounder explains 5% variation of the outcome
        # SNP explains 5% variation of the outcome
        Y <- C*.19 + X*.395 + rnorm(n_obs)
        fit <- lm(Y ~ X + C)
        px[i] <- fit %>% tidy %>% dplyr::filter(term == "X") %>% dplyr::pull(p.value)
        rx[i] <- summary(lm(Y ~ X))$r.squared
        pc[i] <- fit %>% tidy %>% dplyr::filter(term == "C") %>% dplyr::pull(p.value)
        rc[i] <- summary(lm(Y ~ C))$r.squared
    }
    binom.test(sum(px<0.05), n_sim)
    binom.test(sum(pc<0.05), n_sim)
    t.test(rx)
    t.test(rc)
}

check_confounded_estimate <- function(){
    b_crude <- rep(NA, n_sim)
    b_adj <- rep(NA, n_sim)
    for (i in 1:n_sim){
        C <- rnorm(n_obs)
        X <- sapply(C, function(c) {p <- 1 / (1 + exp(-c)); get_simulated_genotypes(p *.255, 1)})
        Y <- C*.19 + X*.395 + rnorm(n_obs)
        b_crude[i] <- lm(Y ~ X) %>% tidy %>% dplyr::filter(term == "X") %>% dplyr::pull(estimate)
        b_adj[i] <- lm(Y ~ X + C) %>% tidy %>%  dplyr::filter(term == "X") %>% dplyr::pull(estimate)
    }
    print(t.test(b_crude))
    print(t.test(b_adj))
}

check_first_stage_adjusted <- function(){
    # check first-stage model is adjusted for covariates
    res_crude <- data.frame()
    res_adj <- data.frame()
    for (i in 1:n_sim){
        C <- rnorm(n_obs)
        X <- sapply(C, function(c) {p <- 1 / (1 + exp(-c)); get_simulated_genotypes(p *.255, 1)})
        Y <- C*.19 + X*.395 + rnorm(n_obs)
        data <- data.frame(
            S = paste0("S", seq(1, n_obs)),
            X,
            Y,
            C,
            stringsAsFactors=F
        )
        res_crude <- rbind(res_crude, run_models(data))
        res_adj <- rbind(res_adj, run_models(data, covar=c("C")))
    }
    results <- rbind(
        t.test(res_adj$BETA_mu.cpp_bp) %>% tidy %>% dplyr::mutate(term="BETA_mu.cpp_bp"),
        t.test(res_crude$BETA_mu.cpp_bp) %>% tidy %>% dplyr::mutate(term="BETA_mu.cpp_bp"),
        t.test(res_adj$BETA_mu.cpp_bf) %>% tidy %>% dplyr::mutate(term="BETA_mu.cpp_bf"),
        t.test(res_crude$BETA_mu.cpp_bf) %>% tidy %>% dplyr::mutate(term="BETA_mu.cpp_bf"),
        t.test(res_adj$BETA_lad.cpp_bf) %>% tidy %>% dplyr::mutate(term="BETA_lad.cpp_bf"),
        t.test(res_crude$BETA_lad.cpp_bf) %>% tidy %>% dplyr::mutate(term="BETA_lad.cpp_bf")
    )
}

check_second_stage_adjusted <- function(covar=NULL){
    # check second-stage model is adjusted for covariates
    results <- data.frame()
    for (i in 1:n_sim){
        C <- rnorm(n_obs)
        Csq <- C^2
        U <- rnorm(n_obs)
        Usq <- U^2
        CU <- C*U
        CUsq <- CU^2
        X <- sapply(C, function(c) {p <- 1 / (1 + exp(-c)); get_simulated_genotypes(p *.9, 1)})
        Y <- C*.28 + X*.79 + C*U*1.12 + rnorm(n_obs)
        data <- data.frame(
            S = paste0("S", seq(1, n_obs)),
            X,
            Y,
            C,
            U,
            CU,
            Csq, Usq, CUsq,
            stringsAsFactors=F
        )
        #res <- run_models(data, covar=covar)
        #res <- data.frame(
        #    rsq.xc=summary(lm(X ~ C))$r.squared,
        #    rsq.yc=summary(lm(Y ~ C))$r.squared,
        #    rsq.yx=summary(lm(Y ~ X))$r.squared,
        #    rsq.ycu=summary(lm(Y ~ C*U))$r.squared
        #)
        res <- vartest(data$Y, data$X, type=2)$test
        res <- vartest(data$Y, data$X, covar = NULL, covar.var = FALSE, type = 1, x.sq = F)
        results <- rbind(results, res)
    }
    return(results)
}

power <- function(results){
    df <- rbind(
        binom.test(sum(results$P.osca_median < .05), n_sim) %>% tidy %>% dplyr::mutate(term="P.osca_median"),
        binom.test(sum(results$P.osca_mean < .05), n_sim) %>% tidy %>% dplyr::mutate(term="P.osca_mean"),
        binom.test(sum(results$P.cpp_bp < .05), n_sim) %>% tidy %>% dplyr::mutate(term="P.cpp_bp"),
        binom.test(sum(results$P.cpp_bf < .05), n_sim) %>% tidy %>% dplyr::mutate(term="P.cpp_bf"),
        binom.test(sum(results$P.r_bf < .05), n_sim) %>% tidy %>% dplyr::mutate(term="P.r_bf"),
        binom.test(sum(results$P.r_bp < .05), n_sim) %>% tidy %>% dplyr::mutate(term="P.r_bp")
    )
}

results <- check_second_stage_adjusted()
results.adj <- check_second_stage_adjusted(covar=c("C", "U", "CU"))
results.adjsq <- check_second_stage_adjusted(covar=c("C", "U", "CU", "Csq", "Usq", "CUsq"))

df <- power(results)
df.adj <- power(results.adj)
df.adjsq <- power(results.adjsq)