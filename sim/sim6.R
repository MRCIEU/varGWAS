library("broom")
library("cqrReg")
library("dplyr")
library("boot")
library("car")
set.seed(23)

n_obs <- 1000
n_sim <- 200

get_residual <- function(x, y, covar=NULL){
    if (!is.null(covar)){
        X <- as.matrix(cbind(x, covar))
    } else {
        X <- as.matrix(data.frame(x))
    }
    # betas
    fit <- qrfit(X=X, y=y, tau=.5, method="mm")
    b <- rbind(fit$b, fit$beta)
    # predicted
    X <- cbind(rep(1, nrow(X)), X)
    fitted <- X %*% b
    # residual
    d <- y - fitted
    return(as.vector(d))
}

model <- function(data, indices){
    # subset vectors
    data <- data[indices,]
    x <- data$x
    y <- data$y
    # absresi
    d <- abs(get_residual(x, y))
    # dummy SNP
    x <- as.factor(x)
    # second-stage model
    fit2 <- lm(d ~ x)
    # extract coef
    b0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull("estimate")
    b1 <- fit2 %>% tidy %>% dplyr::filter(term == "x1") %>% dplyr::pull("estimate")
    b2 <- fit2 %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull("estimate")
    # difference between 0 vs 1 and 0 vs 2
    return(c(
        (2*b0*b1+b1^2)/(2/pi),
        (2*b0*b2+b2^2)/(2/pi)
    ))
}

results <- data.frame()
for (b in seq(2, 2)){
    for (i in 1:n_sim){
        message(paste0("b:", b, " i:", i))
        # SNP
        x <- rbinom(n_obs, 2, .5)
        # modifier
        u <- rnorm(n_obs)
        # outcome
        y <- x*u*b + rnorm(n_obs)
        # estimate variance effects and boostrap SE
        bs <- boot(data.frame(x, y), model, R=1000, stype="i")
        # extract parameters
        e1 <- bs %>% tidy %>% dplyr::pull("statistic") %>% dplyr::nth(1)
        se1 <- bs %>% tidy %>% dplyr::pull("std.error") %>% dplyr::nth(1)
        e2 <- bs %>% tidy %>% dplyr::pull("statistic") %>% dplyr::nth(2)
        se2 <- bs %>% tidy %>% dplyr::pull("std.error") %>% dplyr::nth(1)
        lci1 <- e1 - (1.96 * se1)
        uci1 <- e1 + (1.96 * se1)
        lci2 <- e2 - (1.96 * se2)
        uci2 <- e2 + (1.96 * se2)
        # store results
        results <- rbind(results, data.frame(
            v1=var(y[x==1]) - var(y[x==0]), # true SNP=0 vs SNP=1 variance
            v2=var(y[x==2]) - var(y[x==0]), # true SNP=0 vs SNP=2 variance
            e1, e2, se1, se2, # estimated SNP=0 vs SNP=1 & 2 variance
            lci1, uci1,
            lci2, uci2,
            b
        ))
    }
}

# check for coverage of CI
results %>% dplyr::group_by(b) %>%
    dplyr::summarize(tidy(binom.test(sum(v1 >= lci1 && v1 <= uci1), n_sim))) # count number of times variance of Y is within 95% CI

results %>% dplyr::group_by(b) %>%
    dplyr::summarize(tidy(binom.test(sum(v2 >= lci2 && v2 <= uci2), n_sim))) # count number of times variance of Y is within 95% CI