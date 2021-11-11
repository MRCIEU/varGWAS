library("dplyr")
library("broom")
library("ggplot2")
library("quantreg")
source("funs.R")
set.seed(2445)

n_sim <- 200
n_obs <- 10000
q <- 0.4

results <- data.frame()
for (b in seq(0, 5, 0.5)){
    for (i in 1:n_sim){
        # simulate covariates
        data <- data.frame(
            S = paste0("S", seq(1, n_obs)),
            X = get_simulated_genotypes(q, n_obs),
            stringsAsFactors=F
        )

        # simulate outcome
        data$Y <- rnorm(n_obs, sd=sqrt(1 + data$X * b))
        
        # run models
        res <- run_osca(data, F)

        # estimate variance effect
        fit1 <- lm(Y ~ X, data=data)
        data$d <- resid(fit1)^2
        fit2 <- lm(d ~ X, data=data)
        fit3 <- lm(d ~ 1, data=data)
        p <- anova(fit2, fit3) %>% tidy %>% dplyr::pull(p.value) %>% dplyr::nth(2)
        b0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(estimate)
        s0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(std.error)
        b1 <- fit2 %>% tidy %>% dplyr::filter(term == "X") %>% dplyr::pull(estimate)
        s1 <- fit2 %>% tidy %>% dplyr::filter(term == "X") %>% dplyr::pull(std.error)

        results <- rbind(results, data.frame(b, b0, b1, s0, s1, p, p_osca=res$P.osca_median, b_osca=res$BETA_x.osca_median, s_osca=res$SE_x.osca_median))
    }
}

results %>%
    dplyr::group_by(b) %>%
    dplyr::summarize(b0=mean(b0), b1=mean(b1), b_osca=mean(b_osca))