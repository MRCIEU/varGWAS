# get P value using LAD-BF with dummy method
dummy <- function(x, y, covar=NULL){
    if (!is.null(covar)){
        X <- as.matrix(cbind(x, covar))
    } else {
        X <- as.matrix(data.frame(x))
    }
    # first-stage fit
    fit <- qrfit(X=X, y=y, tau=.5, method="mm")
    b <- rbind(fit$b, fit$beta)
    # predicted
    X <- cbind(rep(1, nrow(X)), X)
    fitted <- X %*% b
    # residual
    d <- y - fitted
    # abs residual
    d <- abs(as.vector(d))
    # dummy SNP
    x <- as.factor(x)
    # second-stage model
    fit2 <- lm(d ~ x)
    fit_null <- lm(d ~ 1)
    p <- anova(fit_null, fit2) %>% tidy %>% dplyr::pull(p.value) %>% dplyr::nth(2)
    return(p)
}

# get P value using LAD-BF with xsq method
xsq <- function(x, y, covar=NULL){
    if (!is.null(covar)){
        X <- as.matrix(cbind(x, covar))
    } else {
        X <- as.matrix(data.frame(x))
    }
    # first-stage fit
    fit <- qrfit(X=X, y=y, tau=.5, method="mm")
    b <- rbind(fit$b, fit$beta)
    # predicted
    X <- cbind(rep(1, nrow(X)), X)
    fitted <- X %*% b
    # residual
    d <- y - fitted
    # abs residual
    d <- abs(as.vector(d))
    # second-stage model
    xsq <- x^2
    fit2 <- lm(d ~ x + xsq)
    fit_null <- lm(d ~ 1)
    p <- anova(fit_null, fit2) %>% tidy %>% dplyr::pull(p.value) %>% dplyr::nth(2)
    return(p)
}

library("dplyr")
library("broom")
library("ggplot2")
source("funs.R")
set.seed(2445)

n_sim <- 200
n_obs <- 5000

results <- data.frame()
for (i in 1:n_sim){
    # simulate covariates
    data <- data.frame(
        S = paste0("S", seq(1, n_obs)),
        X = get_simulated_genotypes(.4, n_obs),
        stringsAsFactors=F
    )

    # simulate outcome
    data$Y <- rnorm(n_obs, sd=sqrt(2 + data$X * 2))
    
    # run models
    res <- run_osca(data)
    res$dp <- dummy(data$X, data$Y)
    res$xp <- xsq(data$X, data$Y)
    oe <- get_osca_effect(res$dp,.4,n_obs,-1)
    res$oe_beta <- oe[1]
    res$oe_se <- oe[2]

    results <- rbind(results, res)
}