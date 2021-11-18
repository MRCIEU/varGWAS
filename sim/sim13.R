library("dplyr")
library("broom")
library("tidyr")
library("cqrReg")
library("jlst")
source("funs.R")
library("data.table")
set.seed(123)

n_sim <- 200
n_obs <- 1000

# get P value using LAD-BF with xsq method
xsq_p <- function(x, y, covar1=NULL, covar2=NULL){
    if (!is.null(covar1)){
        X <- as.matrix(cbind(x, covar1))
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
    if (!is.null(covar2)){
        X <- data.frame(x, xsq)
        X <- cbind(X, covar2)
        fit_null <- lm(d ~ ., data=covar2)
    } else {
        X <- data.frame(x, xsq)
        fit_null <- lm(d ~ 1)
    }
    fit2 <- lm(d ~ ., data=X)
    p <- anova(fit_null, fit2) %>% tidy %>% dplyr::pull(p.value) %>% dplyr::nth(2)
    return(p)
}

results <- data.frame()
for (i in 1:n_sim){
    # simulate variables

    # genotype with MAF of 0.4
    X <- get_simulated_genotypes(0.4, n_obs)
    # modifier
    U <- rnorm(n_obs)
    Usq <- U^2
    # interaction
    XU <- X*U
    XUsq <- XU^2
    
    # X explains % of Y
    # XU explains % of Y
    # outcome
    Y <- X*.4 + X*U*.3775 + rnorm(n_obs)

    data <- data.frame(
        X,
        Y,
        U,
        XU,
        Usq, XUsq,
        stringsAsFactors=F
    )

    res <- data.frame(
        rsq.yx=summary(lm(Y ~ X))$r.squared,
        rsq.yxu=summary(lm(Y ~ X*U))$r.squared,
        p_adj1=xsq_p(data$X, data$Y, covar1 = NULL, covar2 = NULL),
        p_adj2=xsq_p(data$X, data$Y, covar1 = data %>% dplyr::select("U", "XU"), covar2 = NULL)
    )

    results <- rbind(results, res)
}

# test power w/wo adjustment
binom.test(sum(results$p_adj1 < .05), n_sim) %>% tidy
binom.test(sum(results$p_adj2 < .05), n_sim) %>% tidy