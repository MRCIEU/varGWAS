library("dplyr")
library("broom")
library("ggplot2")
library("boot")
source("funs.R")
set.seed(2445)

n_sim <- 30
n_obs <- 1000
q <- 0.4

# function to obtain regression weights
bs <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  result <- bp_x_model(d)
  return(result)
}

bp_x_model <- function(data){
    fit1 <- lm(Y ~ X, data=data)
    data$d <- resid(fit1)^2
    fit2 <- lm(d ~ X, data=data)
    fit3 <- lm(d ~ 1, data=data)
    p <- anova(fit2, fit3) %>% tidy %>% dplyr::pull(p.value) %>% dplyr::nth(2)
    b0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(estimate)
    s0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(std.error)
    b1 <- fit2 %>% tidy %>% dplyr::filter(term == "X") %>% dplyr::pull(estimate)
    s1 <- fit2 %>% tidy %>% dplyr::filter(term == "X") %>% dplyr::pull(std.error)
    return(c(
        b0,s0,b1,s1,p
    ))
}

results <- data.frame()
for (b in seq(0, 8, 2)){
    for (i in 1:n_sim){
        message(paste0("b:", b, " i:", i))
        # simulate covariates
        df <- data.frame(
            S = paste0("S", seq(1, n_obs)),
            X = get_simulated_genotypes(q, n_obs),
            stringsAsFactors=F
        )

        # simulate outcome
        df$Y <- rnorm(n_obs, sd=sqrt(1 + df$X * b))
        #df$Y <- as.numeric(scale(df$Y))
        
        # run models
        fit_main <- bp_x_model(df)
        fit_boot <- boot(data=df, statistic=bs, R=30) %>% tidy
        fit_osca <- run_osca(df, T)
        oe <- get_osca_effect(fit_osca$P.osca_median,q,n_obs,-1)

        results <- rbind(results, data.frame(
            b,
            b0_main=fit_main[1],
            s0_main=fit_main[2],
            b1_main=fit_main[3],
            s1_main=fit_main[4],
            p_main=fit_main[5],
            b0_boot=fit_boot$statistic[1],
            s0_boot=fit_boot$std.error[1],
            bi0_boot=fit_boot$bias[1],
            b1_boot=fit_boot$statistic[3],
            s1_boot=fit_boot$std.error[3],
            bi1_boot=fit_boot$bias[3],
            p_osca=fit_osca$P.osca_median,
            b_osca=fit_osca$BETA_x.osca_median,
            s_osca=fit_osca$SE_x.osca_median,
            b_osca_r=oe[1],
            s_osca_r=oe[2]
        ))
    }
}