library("dplyr")
library("broom")
library("ggplot2")
source("funs.R")
set.seed(2445)

n_sim <- 50
n_obs <- 1000
q <- 0.4

bp_x_dummy_model <- function(data){
    fit1 <- lm(Y ~ X, data=data)
    data$d <- resid(fit1)^2
    data$f <- as.factor(data$X)
    fit2 <- lm(d ~ f, data=data)
    fit3 <- lm(d ~ 1, data=data)
    v1 <- car::deltaMethod(fit2, "1*b1", parameterNames=c("b0", "b1", "b2"), vcov=sandwich::vcovHC(fit2, type = "HC0"))
    v2 <- car::deltaMethod(fit2, "1*b2", parameterNames=c("b0", "b1", "b2"), vcov=sandwich::vcovHC(fit2, type = "HC0"))
    p <- anova(fit2, fit3) %>% tidy %>% dplyr::pull(p.value) %>% dplyr::nth(2)

    res <- data.frame(
        phi_x1=v1$Estimate,
        se_x1=v1$SE,
        phi_x2=v2$Estimate,
        se_x2=v2$SE,
        phi_p=p
    )

    return(res)
}

bp_x_xsq_model <- function(data){
    data$Xsq <- data$X^2
    fit1 <- lm(Y ~ X, data=data)
    data$d <- resid(fit1)^2
    fit2 <- lm(d ~ X + Xsq, data=data)
    fit3 <- lm(d ~ 1, data=data)
    v1 <- car::deltaMethod(fit2, "1*b1 + 1*b2", parameterNames=c("b0", "b1", "b2"), vcov=sandwich::vcovHC(fit2, type = "HC0"))
    v2 <- car::deltaMethod(fit2, "2*b1 + 4*b2", parameterNames=c("b0", "b1", "b2"), vcov=sandwich::vcovHC(fit2, type = "HC0"))
    p <- anova(fit2, fit3) %>% tidy %>% dplyr::pull(p.value) %>% dplyr::nth(2)

    res <- data.frame(
        phi_x1=v1$Estimate,
        se_x1=v1$SE,
        phi_x2=v2$Estimate,
        se_x2=v2$SE,
        phi_p=p
    )

    return(res)
}

w2l <- function(fit, model, v1, v2){
    res <- data.frame(phi=fit$phi_x1, se=fit$se_x1, b=1, v=v1)
    res <- rbind(res, data.frame(phi=fit$phi_x2, se=fit$se_x2, b=2, v=v2))
    res$model <- model
    return(res)
}

models <- function(df, b){
    v1=var(df$Y[df$X==1]) - var(df$Y[df$X==0])
    v2=var(df$Y[df$X==2]) - var(df$Y[df$X==0])
    bp_x_dummy_fit <- w2l(bp_x_dummy_model(df), "X-dummy", v1, v2)
    bp_x_xsq_fit <- w2l(bp_x_xsq_model(df), "X + X^2", v1, v2)
    
    fit <- rbind(bp_x_dummy_fit, bp_x_xsq_fit)
    return(fit)
}

results <- data.frame()
for (b in seq(0, 6, 1)){
    for (i in 1:n_sim){
        message(paste0("b:", b, " i:",i))
        
        # simulate covariates
        df <- data.frame(
            S = paste0("S", seq(1, n_obs)),
            X = get_simulated_genotypes(q, n_obs),
            U = rnorm(n_obs),
            stringsAsFactors=F
        )

        # simulate linear variance effect
        df$Y <- rnorm(n_obs, sd=sqrt(1 + b * df$X))
        res <- models(df, b)
        res$rel <- "Linear"

        results <- rbind(results, res)

        # simulate interaction variance effect
        df$Y <- df$X * df$U * b + rnorm(n_obs)
        res <- models(df, b)
        res$rel <- "Interaction"

        results <- rbind(results, res)
    }
}

# plots
results$lci <- results$phi - (1.96 * results$se)
results$uci <- results$phi + (1.96 * results$se)
results <- results %>% dplyr::mutate(cond = case_when(b == 1 ~ "var(Y|G==1)", b == 2 ~ "var(Y|G==2)"))
ggplot(results, aes(x=v, y=phi)) +
    geom_point() +
    theme_classic() +
    xlab("Variance") +
    ylab("Estimate") +
    facet_wrap(cond ~ model, scales="free")