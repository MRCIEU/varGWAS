library("dplyr")
library("broom")
library("aod")
library("lmtest")
set.seed(123)

n_obs <- 1000
n_sim <- 1000

results <- data.frame()
for (i in 1:n_sim){
    x <- rbinom(n_obs, 2, .5)
    xsq <- x^2
    u <- rnorm(n_obs)
    y <- x*u*.2 + x*.2 + u*.2 + rnorm(n_obs)
    fit1 <- lm(y ~ x)
    d <- resid(fit1)^2
    fit2 <- lm(d ~ x + xsq)
    fit0 <- lm(d~1)
    p_f <- anova(fit0, fit2) %>% tidy %>% dplyr::select(p.value) %>% tidyr::drop_na() %>% as.numeric
    wald <- wald.test(b = coef(fit2), Sigma = vcov(fit2), Terms = 2:3)
    p_wald <- wald$result$chi2[3] %>% as.numeric
    p_lrt <- lrtest(fit2, fit0) %>% tidy %>% dplyr::select(p.value) %>% tidyr::drop_na() %>% as.numeric
    results <- rbind(results, data.frame(p_f, p_wald, p_lrt))
}