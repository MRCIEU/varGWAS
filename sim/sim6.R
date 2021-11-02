library("broom")
library("multcomp")
library("quantreg")
library("dplyr")
set.seed(23)

n_obs <- 1000
n_sim <- 1000

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
        b0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull("estimate")
        b1 <- fit2 %>% tidy %>% dplyr::filter(term == "x1") %>% dplyr::pull("estimate")
        b2 <- fit2 %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull("estimate")
        s0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull("std.error")
        s1 <- fit2 %>% tidy %>% dplyr::filter(term == "x1") %>% dplyr::pull("std.error")
        s2 <- fit2 %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull("std.error")

        results <- rbind(results, data.frame(
            v1=var(y[x==1]),
            v2=var(y[x==2]),
            e1=b0^2/(2/pi) + (2*b0*b1+b1^2)/(2/pi),
            e2=b0^2/(2/pi) + (2*b0*b2+b2^2)/(2/pi),
            se1=s0^2/(2/pi) + (2*s0*s1+s1^2)/(2/pi),
            se2=s0^2/(2/pi) + (2*s0*s2+s2^2)/(2/pi),
            b
        ))
    }
}

# check for coverage
results$l1 <- results$e1 - (1.96 * results$se1)
results$u1 <- results$e1 + (1.96 * results$se1)
results$l2 <- results$e2 - (1.96 * results$se2)
results$u2 <- results$e2 + (1.96 * results$se2)
results$h1 <- (results$v1 >= results$l1 & results$v1 <= results$u1)
results$h2 <- (results$v2 >= results$l2 & results$v2 <= results$u2)
binom.test(sum(results$h1), nrow(results))
binom.test(sum(results$h2), nrow(results))