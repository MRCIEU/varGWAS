library("dplyr")
library("broom")
library("boot")
source("funs.R")
set.seed(134)

n_sim <- 200
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

# function to obtain regression weights
bs <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  result <- bp_x_model(d)
  return(result)
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
for (b in seq(2,2)){
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
        fit_boot <- boot(data=data, statistic=bs, R=75) %>% tidy

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
            b1_x_boot=fit_boot$statistic[3],
            s1_x_boot=fit_boot$std.error[3],
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

# estimate expected variance effect
expected <- results %>% dplyr::group_by(b) %>% dplyr::summarize(mb=mean((v2-v0)*.5), sb=sd((v2-v0)*.5))
results <- merge(results, expected, "b")

# check coverage for X model
results$lci_x <- results$b1_x - (results$s1_x * 1.96)
results$uci_x <- results$b1_x + (results$s1_x * 1.96)
results %>% 
    dplyr::group_by(b) %>%
    dplyr::summarize(binom.test(sum(lci_x <= mb & uci_x >= mb), n()) %>% tidy)
results$lci_x_boot <- results$b1_x - (results$s1_x_boot * 1.96)
results$uci_x_boot <- results$b1_x + (results$s1_x_boot * 1.96)
results %>% 
    dplyr::group_by(b) %>%
    dplyr::summarize(binom.test(sum(lci_x_boot <= mb & uci_x_boot >= mb), n()) %>% tidy)

# check coverage for OSCA
results$b_osca_var <- results$b_osca / (2/pi)
results$s_osca_var <- results$s_osca / (2/pi)
results$lci_b_osca <- results$b_osca_var - (results$s_osca_var * 1.96)
results$uci_b_osca <- results$b_osca_var + (results$s_osca_var * 1.96)
results %>% 
    dplyr::group_by(b) %>%
    dplyr::summarize(binom.test(sum(lci_b_osca <= mb & uci_b_osca >= mb), n()) %>% tidy)