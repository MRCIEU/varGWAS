library("dplyr")
library("broom")
library("multcomp")
library("quantreg")
set.seed(123)

n_obs <- 100000
n_sim <- 30

results <- data.frame()
for (b in seq(.5, 5, .5)){
    for (i in 1:n_sim){
        # simulate interaction
        x <- rbinom(n_obs, 2, .5)
        u <- rnorm(n_obs)
        y <- x*u*b + rnorm(n_obs)

        # estimate variance for genotype groups
        r <- data.frame(x, u, y) %>% group_by(x) %>% summarize(vy=var(y))
        vi0 <- r$vy[1]
        vi1 <- r$vy[2]
        vi2 <- r$vy[3]
        v1 <- vi1 - vi0
        v2 <- vi2 - vi0

        # fit BF
        fit1 <- suppressWarnings(rq(y ~ x, tau=.5))
        dsq <- abs(resid(fit1))
        fit2 <- lm(dsq ~ factor(x))
        
        # get betas
        b0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% pull(estimate)
        b1 <- fit2 %>% tidy %>% dplyr::filter(term == "factor(x)1") %>% pull(estimate)
        b2 <- fit2 %>% tidy %>% dplyr::filter(term == "factor(x)2") %>% pull(estimate)

        B0 <- b0^2/(2/pi)
        B1 <- (2*b0*b1+b1^2)/(2/pi)
        B2 <- (2*b0*b2+b2^2)/(2/pi)
        
        vi0_est <- B0+B1*0+B2*0
        vi1_est <- B0+B1*1+B2*0
        vi2_est <- B0+B1*0+B2*1

        v1_est <- glht(model=fit2, linfct=paste("x*", (2*b0*b1+b1^2)/(2/pi), " == 0")) %>% tidy %>% dplyr::pull(estimate)

        # calculate expected variance in each groups
        results <- rbind(results, data.frame(b, vi0_est, vi1_est, vi2_est, b0, b1, b2, vi0, vi1, vi2))
    }
}