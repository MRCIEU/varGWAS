library("broom")
library("quantreg")
library("dplyr")
library("boot")
library("car")
set.seed(23)

n_obs <- 1000
n_sim <- 200

model <- function(x, y){
    # LAD regression
    fit <- suppressWarnings(rq(y ~ x, tau=0.5))
    # absresi
    d <- abs(resid(fit))
    # dummy SNP
    x <- as.factor(x)
    # second-stage model
    fit2 <- lm(d ~ x)
    # extract coef
    b <- fit2 %>% tidy %>% dplyr::pull("estimate")
    return(b)
}

get_est0 <- function(d, i){
    b <- model(d$x[i], d$y[i])
    return(b[1]^2/(2/pi))
}
get_est1 <- function(d, i){
    b <- model(d$x[i], d$y[i])
    return(b[1]^2/(2/pi) + (2*b[1]*b[2]+b[2]^2)/(2/pi))
}
get_est2 <- function(d, i){
    b <- model(d$x[i], d$y[i])
    return(b[1]^2/(2/pi) + (2*b[1]*b[3]+b[3]^2)/(2/pi))
}

results <- data.frame()
for (b in seq(2)){
    for (i in 1:n_sim){
        # SNP
        x <- rbinom(n_obs, 2, .6)
        # modifier
        u <- rnorm(n_obs)
        # outcome
        y <- x*u*b + rnorm(n_obs)
        # estimate variance for SNP=1 and boostrap SE
        bs1 <- boot(data.frame(x, y), get_est1, R=1000, stype="i") %>% tidy
        # estimate variance for SNP=2 and boostrap SE
        bs2 <- boot(data.frame(x, y), get_est2, R=1000, stype="i") %>% tidy
        # extract parameters
        e1 <- bs1$statistic
        se1 <- bs1$std.error
        e2 <- bs2$statistic
        se2 <- bs2$std.error
        lci1 <- e1 - (1.96 * se1)
        uci1 <- e1 + (1.96 * se1)
        lci2 <- e2 - (1.96 * se2)
        uci2 <- e2 + (1.96 * se2)
        # store results
        results <- rbind(results, data.frame(
            v1=var(y[x==1]), # true variance of SNP=1 group
            v2=var(y[x==2]), # true variance of SNP=2 group
            e1, e2, se1, se2,
            lci1, uci1,
            lci2, uci2,
            b
        ))
    }
}

# check for coverage of CI
results %>% dplyr::group_by(b) %>%
    dplyr::summarize(tidy(binom.test(sum(v1 >= lci1 & v1 <= uci1), n()))) # count number of times variance of Y is within 95% CI

results %>% dplyr::group_by(b) %>%
    dplyr::summarize(tidy(binom.test(sum(v2 >= lci2 & v2 <= uci2), n()))) # count number of times variance of Y is within 95% CI