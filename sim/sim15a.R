library("broom")
source("funs.R")
set.seed(123)

# MAF
n_obs <- 10000
q <- 0.25

# effect of var(Y|X) under a range of conditions

# main effect of X has no effect on var(Y)
results <- data.frame()
for (b in seq(0, 50, 0.25)){
    x <- get_simulated_genotypes(q, n_obs)
    y <- x*b + rnorm(n_obs)
    v1 <- var(y[x==1]) - var(y[x==0])
    v2 <- var(y[x==2]) - var(y[x==0])
    results <- rbind(results, data.frame(b, v1, v2))
}
plot(results$v1 ~ results$b)
lm(v1 ~ b, data=results) %>% confint
plot(results$v2 ~ results$b)
lm(v2 ~ b, data=results) %>% confint

# interaction effect on Y
# X^2 correlates with var(Y)
# model = var(Y|X==1) X^2 * 1
# model = var(Y|X==2) X^2 * 4
results <- data.frame()
for (b in seq(0, 50, 0.25)){
    x <- get_simulated_genotypes(q, n_obs)
    u <- rnorm(n_obs)
    y <- x*u*b + rnorm(n_obs)
    v1 <- var(y[x==1]) - var(y[x==0])
    v2 <- var(y[x==2]) - var(y[x==0])
    results <- rbind(results, data.frame(b, v1, v2))
}
results$b2 <- results$b^2
results$b2_4 <- results$b2 * 4
plot(results$v1 ~ results$b2)
plot(results$v2 ~ results$b2_4)
lm(v1 ~ b2, data=results) %>% summary
lm(v2 ~ b2_4, data=results) %>% summary


# interaction effect on Y with regression estimate of X
results <- data.frame()
for (b in seq(0, 50, 0.25)){
    x <- get_simulated_genotypes(q, n_obs)
    u <- rnorm(n_obs)
    y <- x*u*b + rnorm(n_obs)
    fit <- lm(y ~ x)
    x2 <- x^2
    d <- resid(fit)^2
    fit <- lm(d ~ x2)
    est <- fit %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull(estimate)
    v <- var(y[x==2]) - var(y[x==0])
    results <- rbind(results, data.frame(b, v, est))
}
plot(v ~ est ,data=results)
lm(v ~ est, data=results) %>% summary


# interaction effect on Y + main effect of X
# X^2 correlates with var(Y)
# model = var(Y|X==1) X^2 * 1
# model = var(Y|X==2) X^2 * 4
results <- data.frame()
for (b in seq(0, 50, 0.25)){
    x <- get_simulated_genotypes(q, n_obs)
    u <- rnorm(n_obs)
    y <- x*b + x*u*b + rnorm(n_obs)
    v1 <- var(y[x==1]) - var(y[x==0])
    v2 <- var(y[x==2]) - var(y[x==0])
    results <- rbind(results, data.frame(b, v1, v2))
}
results$b2 <- results$b^2
results$b2_4 <- results$b2 * 4
plot(results$v1 ~ results$b2)
plot(results$v2 ~ results$b2_4)
lm(v1 ~ b2, data=results) %>% summary
lm(v2 ~ b2_4, data=results) %>% summary

# interaction effect on Y + main effect of U
# modifier has effect on var(Y)
results <- data.frame()
for (b in seq(0, 50, 0.25)){
    x <- get_simulated_genotypes(q, n_obs)
    u <- rnorm(n_obs)
    y <- u*b + x*u*b + rnorm(n_obs)
    v1 <- var(y[x==1]) - var(y[x==0])
    v2 <- var(y[x==2]) - var(y[x==0])
    results <- rbind(results, data.frame(b, v1, v2))
}
results$b2 <- results$b^2
results$b2_4 <- results$b2 * 4
plot(results$v1 ~ results$b2)
plot(results$v2 ~ results$b2_4)
lm(v1 ~ b2, data=results) %>% summary
lm(v2 ~ b2_4, data=results) %>% summary

# interaction effect on Y
# regression estimate of X on var(Y)
results <- data.frame()
for (b in seq(0, 50, 0.25)){
    x <- get_simulated_genotypes(q, n_obs)
    u <- rnorm(n_obs)
    y <- x*u*b + rnorm(n_obs)
    v1 <- var(y[x==1]) - var(y[x==0])
    v2 <- var(y[x==2]) - var(y[x==0])
    fit <- lm(y ~ x)
    d <- resid(fit)^2
    x2 <- x^2
    fit <- lm(d ~ x2)
    v1_est <- fit %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull(estimate) * 1
    v2_est <- fit %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull(estimate) * 4
    results <- rbind(results, data.frame(b, v1, v2, v1_est, v2_est))
}
plot(results$v1 ~ results$v1_est)
plot(results$v2 ~ results$v2_est)

# interaction effect on Y + main effect of X
# regression estimate of X on var(Y)
results <- data.frame()
for (b in seq(0, 50, 0.25)){
    x <- get_simulated_genotypes(q, n_obs)
    u <- rnorm(n_obs)
    y <- x*b + x*u*b + rnorm(n_obs)
    v1 <- var(y[x==1]) - var(y[x==0])
    v2 <- var(y[x==2]) - var(y[x==0])
    fit <- lm(y ~ x)
    d <- resid(fit)^2
    x2 <- x^2
    fit <- lm(d ~ x2)
    v1_est <- fit %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull(estimate) * 1
    v2_est <- fit %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull(estimate) * 4
    results <- rbind(results, data.frame(b, v1, v2, v1_est, v2_est))
}
plot(results$v1 ~ results$v1_est)
plot(results$v2 ~ results$v2_est)

# interaction effect on Y + main effect of U
# regression estimate of X on var(Y)
# need to include X + X^2 in the second-stage model otherwise if U has a main effect X^2 only will be wrong
results <- data.frame()
for (b in seq(0, 50, 0.25)){
    x <- get_simulated_genotypes(q, n_obs)
    u <- rnorm(n_obs)
    y <- u*b + x*u*b + rnorm(n_obs)
    v1 <- var(y[x==1]) - var(y[x==0])
    v2 <- var(y[x==2]) - var(y[x==0])
    fit <- lm(y ~ x)
    d <- resid(fit)^2
    x2 <- x^2
    fit <- lm(d ~ x + x2)
    v1_est <- fit %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::pull(estimate) * 1 +
        fit %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull(estimate) * 1
    v2_est <- fit %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::pull(estimate) * 2 +
        fit %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull(estimate) * 4
    results <- rbind(results, data.frame(b, v1, v2, v1_est, v2_est))
}
plot(results$v1 ~ results$v1_est)
plot(results$v2 ~ results$v2_est)