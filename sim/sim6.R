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
        x <- rbinom(n_obs, 2, .5)
        u <- rnorm(n_obs)
        y <- x*u*b + rnorm(n_obs)

        d <- data.frame(x, u, y)
        r <- d %>% group_by(x) %>% summarize(my=mean(y), vy=var(y), mady=mad(y))
        m1 <- r$my[2] - r$my[1]
        m2 <- r$my[3] - r$my[1]
        v1 <- r$vy[2] - r$vy[1]
        v2 <- r$vy[3] - r$vy[1]
        s1 <- sqrt(v1)
        s2 <- sqrt(v2)
        mad1 <- r$mady[2] - r$mady[1]
        mad2 <- r$mady[3] - r$mady[1]

        # BP
        fit1 <- lm(y ~ x)
        dsq <- resid(fit1)^2
        xsq <- x^2
        fit2 <- lm(dsq ~ x + xsq)
        
        varbeta1 <- glht(model=fit2, linfct=paste("x*1 + xsq*1 == 0")) %>% tidy %>% dplyr::pull(estimate)
        varbeta2 <- glht(model=fit2, linfct=paste("x*2 + xsq*4 == 0")) %>% tidy %>% dplyr::pull(estimate)

        # BF
        fit1 <- suppressWarnings(rq(y ~ x, tau=.5))
        dsq <- abs(resid(fit1))
        xsq <- x^2
        fit2 <- lm(dsq ~ x + xsq)
        
        madbeta1 <- glht(model=fit2, linfct=paste("x*", 1/sqrt(2/pi), " + xsq*1 == 0")) %>% tidy %>% dplyr::pull(estimate)
        madbeta2 <- glht(model=fit2, linfct=paste("x*", 2*1/sqrt(2/pi), " + xsq*", 4*1/sqrt(2/pi), " == 0")) %>% tidy %>% dplyr::pull(estimate)
        sigmabeta1 <- glht(model=fit2, linfct=paste("x*1 + xsq*1 == 0")) %>% tidy %>% dplyr::pull(estimate)
        sigmabeta2 <- glht(model=fit2, linfct=paste("x*2 + xsq*4 == 0")) %>% tidy %>% dplyr::pull(estimate)

        int <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% pull(estimate)
        xb <- fit2 %>% tidy %>% dplyr::filter(term == "x") %>% pull(estimate)
        xsqb <- fit2 %>% tidy %>% dplyr::filter(term == "xsq") %>% pull(estimate)

        results <- rbind(results, data.frame(b, m1, m2, v1, v2, s1, s2, mad1, mad2, varbeta1, varbeta2, madbeta1, madbeta2, sigmabeta1, sigmabeta2, int, xb, xsqb))
    }
}

# BP
v1 <- results %>% dplyr::group_by(b) %>% dplyr::summarize(t.test(v1) %>% tidy)
v2 <- results %>% dplyr::group_by(b) %>% dplyr::summarize(t.test(v2) %>% tidy)
varbeta1 <- results %>% dplyr::group_by(b) %>% dplyr::summarize(t.test(varbeta1) %>% tidy)
varbeta2 <- results %>% dplyr::group_by(b) %>% dplyr::summarize(t.test(varbeta2) %>% tidy)

# BF
mad1 <- results %>% dplyr::group_by(b) %>% dplyr::summarize(t.test(mad1) %>% tidy)
mad2 <- results %>% dplyr::group_by(b) %>% dplyr::summarize(t.test(mad2) %>% tidy)
madbeta1 <- results %>% dplyr::group_by(b) %>% dplyr::summarize(t.test(madbeta1) %>% tidy)
madbeta2 <- results %>% dplyr::group_by(b) %>% dplyr::summarize(t.test(madbeta2) %>% tidy)
sigmabeta1 <- results %>% dplyr::group_by(b) %>% dplyr::summarize(t.test(sigmabeta1) %>% tidy)
sigmabeta2 <- results %>% dplyr::group_by(b) %>% dplyr::summarize(t.test(sigmabeta2) %>% tidy)