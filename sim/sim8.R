library("dplyr")
library("broom")
source("funs.R")
set.seed(134)

n_sim <- 30
n_obs <- 1000
q <- 0.4

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
    p1 <- fit2 %>% tidy %>% dplyr::filter(term == "X") %>% dplyr::pull(p.value)
    return(c(
        b0,s0,b1,s1,p1
    ))
}

bp_xabs_model <- function(data){
    fit1 <- suppressWarnings(rq(Y ~ X, data=data))
    data$d <- abs(resid(fit1))
    fit2 <- lm(d ~ X, data=data)
    b0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(estimate)
    s0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(std.error)
    b1 <- fit2 %>% tidy %>% dplyr::filter(term == "X") %>% dplyr::pull(estimate)
    s1 <- fit2 %>% tidy %>% dplyr::filter(term == "X") %>% dplyr::pull(std.error)
    p1 <- fit2 %>% tidy %>% dplyr::filter(term == "X") %>% dplyr::pull(p.value)
    return(c(
        b0,s0,b1,s1,p1
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
for (b in seq(0, 8, 2)){
    for (i in 1:n_sim){
        message(paste0("b:", b, " i:", i))
        # simulate covariates
        data <- data.frame(
            S = paste0("S", seq(1, n_obs)),
            X = get_simulated_genotypes(q, n_obs),
            U = rnorm(n_obs),
            stringsAsFactors=F
        )

        # simulate outcome
        data$Y <- rnorm(n_obs, sd=sqrt(1 + data$X * b))
        #data$Y <- data$X * data$U * b + rnorm(n_obs)

        # standardise trait
        data$Y <- scale(data$Y)
        
        # fit models
        fit_osca <- run_osca(data, T)
        fit_xd <- bp_x_dummy_model(data)
        fit_xsq <- bp_xsq_model(data)
        fit_x <- bp_x_model(data)
        fit_xabs <- bp_xabs_model(data)

        results <- rbind(results, data.frame(
            b,
            b0_xsq=fit_xsq[1],
            s0_xsq=fit_xsq[2],
            b1_xsq=fit_xsq[3],
            s1_xsq=fit_xsq[4],
            b2_xsq=fit_xsq[5],
            s2_xsq=fit_xsq[6],
            p_xsq=fit_xsq[7],
            b0_xd=fit_xd[1],
            s0_xd=fit_xd[2],
            b1_xd=fit_xd[3],
            s1_xd=fit_xd[4],
            b2_xd=fit_xd[5],
            s2_xd=fit_xd[6],
            p_xd=fit_xd[7],
            b0_x=fit_x[1],
            s0_x=fit_x[2],
            b1_x=fit_x[3],
            s1_x=fit_x[4],
            p_x=fit_x[5],
            b0_xabs=fit_xabs[1],
            s0_xabs=fit_xabs[2],
            b1_xabs=fit_xabs[3],
            s1_xabs=fit_xabs[4],
            p_xabs=fit_xabs[5],
            v0=var(data$Y[data$X==0]),
            v1=var(data$Y[data$X==1]),
            v2=var(data$Y[data$X==2]),
            mad0=mad(data$Y[data$X==0]),
            mad1=mad(data$Y[data$X==1]),
            mad2=mad(data$Y[data$X==2]),
            b_osca=fit_osca$BETA_x.osca_median,
            s_osca=fit_osca$SE_x.osca_median,
            p_osca=fit_osca$P.osca_median,
            n_obs, q
        ))

    }
}