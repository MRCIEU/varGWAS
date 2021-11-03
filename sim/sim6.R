library("broom")
library("quantreg")
library("dplyr")
library("car")
set.seed(23)

n_obs <- 1000
n_sim <- 200

results <- data.frame()
for (b in seq(0, 6, .5)){
    for (i in 1:n_sim){
        x <- rbinom(n_obs, 2, .6)
        u <- rnorm(n_obs)
        y <- x*u*b + rnorm(n_obs)

        fit <- suppressWarnings(rq(y ~ x, tau=0.5))
        d <- abs(resid(fit))
        x <- as.factor(x)
        fit2 <- lm(d ~ x)
        d1 <- deltaMethod(fit2, "(Intercept)^2/(2/pi) + (2*(Intercept)*x1+x1^2)/(2/pi)", vcov=hccm)
        d2 <- deltaMethod(fit2, "(Intercept)^2/(2/pi) + (2*(Intercept)*x2+x2^2)/(2/pi)", vcov=hccm)
        e1 <- d1 %>% dplyr::pull("Estimate")
        e2 <- d2 %>% dplyr::pull("Estimate")
        se1 <- d1 %>% dplyr::pull("SE")
        se2 <- d2 %>% dplyr::pull("SE")
        lci1 <- d1 %>% dplyr::pull("2.5 %")
        uci1 <- d1 %>% dplyr::pull("97.5 %")
        lci2 <- d2 %>% dplyr::pull("2.5 %")
        uci2 <- d2 %>% dplyr::pull("97.5 %")

        results <- rbind(results, data.frame(
            v1=var(y[x==1]),
            v2=var(y[x==2]),
            e1, e2, se1, se2,
            lci1, uci1,
            lci2, uci2,
            b
        ))
    }
}

# check for coverage of CI
results %>% dplyr::group_by(b) %>%
    dplyr::summarize(tidy(binom.test(sum(v1 >= lci1 & v1 <= uci1), n())))

results %>% dplyr::group_by(b) %>%
    dplyr::summarize(tidy(binom.test(sum(v2 >= lci2 & v2 <= uci2), n())))