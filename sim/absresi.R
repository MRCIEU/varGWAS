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
        r <- data.frame(x, u, y) %>% group_by(x) %>% summarize(vy=var(y), mad=mad(y))
        vi0 <- r$vy[1]
        vi1 <- r$vy[2]
        vi2 <- r$vy[3]
        mi0 <- r$mad[1]
        mi1 <- r$mad[2]
        mi2 <- r$mad[3]
        v1 <- vi1 - vi0
        v2 <- vi2 - vi0

        # fit BF
        fit1 <- suppressWarnings(rq(y ~ x, tau=.5))
        dsq <- abs(resid(fit1))
        xsq <- x^2
        fit2 <- lm(dsq ~ x + xsq)
        
        # get betas
        b0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% pull(estimate)
        b1 <- fit2 %>% tidy %>% dplyr::filter(term == "x") %>% pull(estimate)
        b2 <- fit2 %>% tidy %>% dplyr::filter(term == "xsq") %>% pull(estimate)

        # calculate expected variance in each group
        # MAD = SD * sqrt(2/pi)
        # var = MAD^2 / (2/pi)
        B0 <- b0^2/(2/pi)
        B1 <- 2*b0*b1/(2/pi)
        B2 <- b1^2/(2/pi)
        varbetai0_exp <- B0+B1*0+B2*0^2
        varbetai1_exp <- B0+B1*1+B2*1^2
        varbetai2_exp <- B0+B1*2+B2*4^2
        varbeta1_exp <- B1*1+B2*1^2
        varbeta2_exp <- B1*2+B2*4^2

        results <- rbind(results, data.frame(b, v1, v2, b0, b1, b2, varbeta1_exp, varbeta2_exp, varbetai0_exp, varbetai1_exp, varbetai2_exp, vi0, vi1, vi2, mi0, mi1, mi2))
    }
}

plot(results$varbeta1_exp, results$v1)
plot(results$varbeta2_exp, results$v2)

results %>% 
    group_by(b) %>% 
    summarize(
        m=mean(
            (
                b0 * 1 / sqrt(2/pi) +
                b1 * 2 / sqrt(2/pi) +
                b2 * 4 / sqrt(2/pi)
            )^2
        ), v=mean(vi2)) %>%
        with(plot(m,v))

tx <- results %>% 
    group_by(b) %>% 
    summarize(
        m=mean(
            (
                b0^2 / (2/pi) +
                2*b0*b1 / (2/pi) * 1 +
                2*b0*b2 / (2/pi) * 1 +
                2*b1*b2 / (2/pi) * 1
                b1^2 / (2/pi) * 1 +
                b2^2 / (2/pi) * 1
            )
        ), v=mean(vi1)
    ) %>% 
    do(tidy(lm(v~m, .)))
    with(plot(m,v))

tx <- results %>% 
    group_by(b) %>% 
    summarize(
        m=mean(
            (
                b0^2 / (2/pi) +
                ((b1^2 / (2/pi)) + (2*b0*b1 / (2/pi)) + (2*b1*b2 / (2/pi))) * 2 +
                ((b2^2 / (2/pi)) + (2*b0*b2 / (2/pi)) + (2*b1*b2 / (2/pi))) * 4
            )
        ), v=mean(vi2)
    ) %>% 
    with(plot(m,v))