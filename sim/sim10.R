library("dplyr")
library("broom")
library("ggplot2")
library("boot")
source("funs.R")
set.seed(2445)

n_sim <- 50
n_obs <- 1000
q <- 0.4

# function to obtain regression weights
bs <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  result <- model(d)
  return(result)
}

model <- function(data){
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
for (b in c(0, 5)){
    for (i in 1:n_sim){
        message(paste0("b:", b, " i:", i))
        # simulate covariates
        data <- data.frame(
            S = paste0("S", seq(1, n_obs)),
            X = get_simulated_genotypes(q, n_obs),
            stringsAsFactors=F
        )

        # simulate outcome
        data$Y <- rnorm(n_obs, sd=sqrt(1 + data$X * b))
        
        # run models
        fit_main <- model(data)
        fit_boot <- boot(data=data, statistic=bs, R=50) %>% tidy
        fit_osca <- run_osca(data, T)

        results <- rbind(results, data.frame(
            b,
            b0_main=fit_main[1],
            s0_main=fit_main[3],
            b1_main=fit_main[2],
            s1_main=fit_main[4],
            p_main=fit_main[5],
            b0_boot=fit_boot$statistic[1],
            s0_boot=fit_boot$std.error[1],
            b1_boot=fit_boot$statistic[2],
            s1_boot=fit_boot$std.error[2],
            p_osca=res$P.osca_median,
            b_osca=res$BETA_x.osca_median,
            s_osca=res$SE_x.osca_median
        ))
    }
}