library("dplyr")
library("broom")
library("ggplot2")
library("multcomp")
library("boot")
source("funs.R")
set.seed(2445)

n_sim <- 50
n_obs <- 1000
q <- 0.4

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

bp_x_dummy_model <- function(data){
    fit1 <- lm(Y ~ X, data=data)
    data$d <- resid(fit1)^2
    data$X <- as.factor(data$X)
    fit2 <- lm(d ~ X, data=data)
    fit3 <- lm(d ~ 1, data=data)
    p <- anova(fit2, fit3) %>% tidy %>% dplyr::pull(p.value) %>% dplyr::nth(2)
    b0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(estimate)
    s0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(std.error)
    b1 <- fit2 %>% tidy %>% dplyr::filter(term == "X1") %>% dplyr::pull(estimate)
    s1 <- fit2 %>% tidy %>% dplyr::filter(term == "X1") %>% dplyr::pull(std.error)
    b2 <- fit2 %>% tidy %>% dplyr::filter(term == "X2") %>% dplyr::pull(estimate)
    s2 <- fit2 %>% tidy %>% dplyr::filter(term == "X2") %>% dplyr::pull(std.error)
    return(c(
        b0,s0,b1,s1,b2,s2,p
    ))
}

bp_xsq_model <- function(data){
    data$xsq <- data$X^2
    fit1 <- lm(Y ~ X, data=data)
    data$d <- resid(fit1)^2
    fit2 <- lm(d ~ X + xsq, data=data)
    fit3 <- lm(d ~ 1, data=data)
    p <- anova(fit2, fit3) %>% tidy %>% dplyr::pull(p.value) %>% dplyr::nth(2)
    b0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(estimate)
    s0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(std.error)
    b1 <- fit2 %>% tidy %>% dplyr::filter(term == "X") %>% dplyr::pull(estimate)
    s1 <- fit2 %>% tidy %>% dplyr::filter(term == "X") %>% dplyr::pull(std.error)
    b2 <- fit2 %>% tidy %>% dplyr::filter(term == "xsq") %>% dplyr::pull(estimate)
    s2 <- fit2 %>% tidy %>% dplyr::filter(term == "xsq") %>% dplyr::pull(std.error)
    return(c(
        b0,s0,b1,s1,b2,s2,p
    ))
}

results <- data.frame()
for (b in seq(0, 6, 1)){
    for (i in 1:n_sim){
        message(paste0("b:", b, " i:", i))
        
        # simulate covariates
        df <- data.frame(
            S = paste0("S", seq(1, n_obs)),
            X = get_simulated_genotypes(q, n_obs),
            U = rnorm(n_obs),
            stringsAsFactors=F
        )

        # simulate outcome
        df$Y <- df$X * df$U * b + rnorm(n_obs)
        #df$Y <- rnorm(n_obs, sd=sqrt(1 + b*df$X))
        
        # run models
        fit_main <- bp_x_dummy_model(df)

        results <- rbind(results, data.frame(
            b,
            b0_main=fit_main[1],
            s0_main=fit_main[2],
            b1_main=fit_main[3],
            s1_main=fit_main[4],
            b2_main=fit_main[5],
            s2_main=fit_main[6],
            p_main=fit_main[7],
            v0=var(df$Y[df$X==0]),
            v1=var(df$Y[df$X==1]),
            v2=var(df$Y[df$X==2])
        ))
    }
}

# include X only in the second-stage model = var(Y|G) is not linear when an interaction is present but is linear if Y is sampled from sigma*b
#results %>% dplyr::group_by(b) %>% dplyr::summarize(v0m=mean(v0), v1m=mean(v1), v2m=mean(v2))

# include X + xsq and you will get the correct value for variance if
# varince g1=x*1+xsq*1
# varince g2=x*2+xsq*4
#r <- results %>% dplyr::group_by(b) %>% dplyr::summarize(b1m=mean(b1_main), b2m=mean(b2_main), v1m=mean(v1), v2m=mean(v2))

# dummy method also gives correct value for variance
r <- results %>% dplyr::group_by(b) %>% dplyr::summarize(b1m=mean(b1_main), b2m=mean(b2_main), v1m=mean(v1), v2m=mean(v2))