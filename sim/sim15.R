library("broom")
library("ggplot2")
library("dplyr")
source("funs.R")
set.seed(123)

n_sim <- 200
n_obs <- 1000
q <- 0.25 # MAF
b <- 0.145 # main effect size to have 80% power

# interaction effect on change in var(Y)
results <- data.frame()
for (i in 1:n_sim){
    for (phi in seq(0, 6, 0.25)){
        theta <- phi * b
        x <- get_simulated_genotypes(q, n_obs)
        u <- rnorm(n_obs)
        y <- x*b + x*u*theta + u*0 + rnorm(n_obs)
        results <- rbind(results, data.frame(
            phi, var=var(y) - 1, u="None"
        ))
        y <- x*b + x*u*theta + u*b + rnorm(n_obs)
        results <- rbind(results, data.frame(
            phi, var=var(y) - 1, u="Detectable with 80% power"
        ))
    }
}

summary <- results %>%
    dplyr::group_by(phi, u) %>% 
    dplyr::summarise(t.test(var) %>% tidy)

pdf("Interaction_variance.pdf")
ggplot(summary, aes(x=phi, y=estimate, ymin=conf.low, ymax=conf.high, group=u, shape=u)) +
    geom_point() + geom_errorbar(width=.05) + theme_classic() +
    labs(y="Change in total variance (95% CI)",x="Interaction effect size relative to main effect",shape="Modifier main effect")
dev.off()

# genotype effect of interaction on var(Y)
results <- data.frame()
for (i in 1:n_sim) {
    theta <- 6 * b
    x <- get_simulated_genotypes(q, n_obs)
    u <- rnorm(n_obs)
    y <- x*u*theta + rnorm(n_obs)
    res <- data.frame(
        v0=var(y[x==0]),
        v1=var(y[x==1]),
        v2=var(y[x==2]),
        model="XU"
    )
    results <- rbind(res, results)
    y <- x*b + x*u*theta + rnorm(n_obs)
    res <- data.frame(
        v0=var(y[x==0]),
        v1=var(y[x==1]),
        v2=var(y[x==2]),
        model="X + XU"
    )
    results <- rbind(res, results)
    y <- u*b + x*u*theta + rnorm(n_obs)
    res <- data.frame(
        v0=var(y[x==0]),
        v1=var(y[x==1]),
        v2=var(y[x==2]),
        model="U + XU"
    )
    results <- rbind(res, results)
}

estimates0 <- results %>% 
    dplyr::group_by(model) %>% 
    dplyr::summarise(t.test(v0) %>% tidy %>% dplyr::mutate(genotype="0"))
estimates1 <- results %>% 
    dplyr::group_by(model) %>% 
    dplyr::summarise(t.test(v1) %>% tidy %>% dplyr::mutate(genotype="1"))
estimates2 <- results %>% 
    dplyr::group_by(model) %>% 
    dplyr::summarise(t.test(v2) %>% tidy %>% dplyr::mutate(genotype="2"))
estimates <- rbind(
    estimates0,estimates1,estimates2
)
estimates$model <- factor(estimates$model, levels=c("XU", "X + XU", "U + XU"))

pdf("Genotype_variance.pdf")
ggplot(estimates, aes(x=genotype, y=estimate, ymin=conf.low, ymax=conf.high, group=model, shape=model)) +
    geom_point(position=position_dodge(width=0.3)) + geom_errorbar(width=.05) + theme_classic() + 
    labs(y="Variance (95% CI)",x="Genotype",shape="Model")
dev.off()

# define quadratic variance model
results <- data.frame()
for (i in 1:n_sim){
    for (phi in seq(0, 6, 0.25)){
        theta <- phi * b
        x <- get_simulated_genotypes(q, n_obs)
        x2 <- x^2
        u <- rnorm(n_obs)
        y <- x*u*theta + rnorm(n_obs)
        fit <- lm(y ~ x)
        d <- resid(fit)
        d2 <- d^2
        res <- data.frame(
            var0=var(y[x==0]),
            var1=var(y[x==1]),
            var2=var(y[x==2])
        )
        fit <- lm(d2 ~ x2)
        b0 <- fit %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(estimate)
        b1 <- fit %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull(estimate)
        res$e0 <- b0
        res$e1 <- b0 + b1 * 1
        res$e2 <- b0 + b1 * 4
        res$phi <- phi

        results <- rbind(results, res)
    }
}

e0 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e0) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==0)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
e1 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e1) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==1)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
e2 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e2) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==2)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
v0 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var0))
v1 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var1))
v2 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var2))
e0 <- cbind(e0, v0 %>% dplyr::select(var))
e1 <- cbind(e1, v1 %>% dplyr::select(var))
e2 <- cbind(e2, v2 %>% dplyr::select(var))
e <- rbind(e0,e1,e2)

pdf("Quadratic_model.pdf")
ggplot(e, aes(x=var, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() + geom_errorbar(width=.05) + theme_classic() + geom_abline(intercept = 0, slope = 1, linetype="dashed", color="grey") + facet_wrap(~genotype, scales="free") +
    labs(y="Estimated variance (95% CI)",x="Outcome variance conditional on genotype")
dev.off()

# modifier main effect with quadratic model
results <- data.frame()
for (i in 1:n_sim){
    for (phi in seq(0, 6, 0.25)){
        theta <- phi * b
        x <- get_simulated_genotypes(q, n_obs)
        x2 <- x^2
        u <- rnorm(n_obs)
        y <- u*theta + x*u*theta + rnorm(n_obs)
        fit <- lm(y ~ x)
        d <- resid(fit)
        d2 <- d^2
        res <- data.frame(
            var0=var(y[x==0]),
            var1=var(y[x==1]),
            var2=var(y[x==2])
        )
        fit <- lm(d2 ~ x2)
        b0 <- fit %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(estimate)
        b1 <- fit %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull(estimate)
        res$e0 <- b0
        res$e1 <- b0 + b1 * 1
        res$e2 <- b0 + b1 * 4
        res$phi <- phi

        results <- rbind(results, res)
    }
}

e0 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e0) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==0)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
e1 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e1) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==1)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
e2 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e2) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==2)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
v0 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var0))
v1 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var1))
v2 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var2))
e0 <- cbind(e0, v0 %>% dplyr::select(var))
e1 <- cbind(e1, v1 %>% dplyr::select(var))
e2 <- cbind(e2, v2 %>% dplyr::select(var))
e <- rbind(e0,e1,e2)

pdf("Quadratic_model_modifier_main.pdf")
ggplot(e, aes(x=var, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() + geom_errorbar(width=.05) + theme_classic() + geom_abline(intercept = 0, slope = 1, linetype="dashed", color="grey") + facet_wrap(~genotype, scales="free") +
    labs(y="Estimated variance (95% CI)",x="Outcome variance conditional on genotype")
dev.off()

# modifier main effect with X + X^2 model
results <- data.frame()
for (i in 1:n_sim){
    for (phi in seq(0, 6, 0.25)){
        theta <- phi * b
        x <- get_simulated_genotypes(q, n_obs)
        x2 <- x^2
        u <- rnorm(n_obs)
        y <- u*theta + x*u*theta + rnorm(n_obs)
        fit <- lm(y ~ x)
        d <- resid(fit)
        d2 <- d^2
        res <- data.frame(
            var0=var(y[x==0]),
            var1=var(y[x==1]),
            var2=var(y[x==2])
        )
        fit <- lm(d2 ~ x + x2)
        b0 <- fit %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(estimate)
        b1 <- fit %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::pull(estimate)
        b2 <- fit %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull(estimate)
        res$e0 <- b0
        res$e1 <- b0 + b1 * 1 + b2 * 1
        res$e2 <- b0 + b1 * 2 + b2 * 4
        res$phi <- phi

        results <- rbind(results, res)
    }
}

e0 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e0) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==0)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
e1 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e1) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==1)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
e2 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e2) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==2)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
v0 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var0))
v1 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var1))
v2 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var2))
e0 <- cbind(e0, v0 %>% dplyr::select(var))
e1 <- cbind(e1, v1 %>% dplyr::select(var))
e2 <- cbind(e2, v2 %>% dplyr::select(var))
e <- rbind(e0,e1,e2)

pdf("Quadratic_model_modifier_main.pdf")
ggplot(e, aes(x=var, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() + geom_errorbar(width=.05) + theme_classic() + geom_abline(intercept = 0, slope = 1, linetype="dashed", color="grey") + facet_wrap(~genotype, scales="free") +
    labs(y="Estimated variance (95% CI)",x="Outcome variance conditional on genotype")
dev.off()

# modifier main effect with dummy method
results <- data.frame()
for (i in 1:n_sim){
    for (phi in seq(0, 6, 0.25)){
        theta <- phi * b
        x <- get_simulated_genotypes(q, n_obs)
        u <- rnorm(n_obs)
        y <- u*theta + x*u*theta + rnorm(n_obs)
        fit <- lm(y ~ x)
        d <- resid(fit)
        d2 <- d^2
        res <- data.frame(
            var0=var(y[x==0]),
            var1=var(y[x==1]),
            var2=var(y[x==2])
        )
        fit <- lm(d2 ~ as.factor(x))
        b0 <- fit %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(estimate)
        b1 <- fit %>% tidy %>% dplyr::filter(term == "as.factor(x)1") %>% dplyr::pull(estimate)
        b2 <- fit %>% tidy %>% dplyr::filter(term == "as.factor(x)2") %>% dplyr::pull(estimate)
        res$e0 <- b0
        res$e1 <- b0 + b1
        res$e2 <- b0 + b2
        res$phi <- phi

        results <- rbind(results, res)
    }
}

e0 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e0) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==0)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
e1 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e1) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==1)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
e2 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e2) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==2)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
v0 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var0))
v1 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var1))
v2 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var2))
e0 <- cbind(e0, v0 %>% dplyr::select(var))
e1 <- cbind(e1, v1 %>% dplyr::select(var))
e2 <- cbind(e2, v2 %>% dplyr::select(var))
e <- rbind(e0,e1,e2)

pdf("Dummy_model_modifier_main.pdf")
ggplot(e, aes(x=var, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() + geom_errorbar(width=.05) + theme_classic() + geom_abline(intercept = 0, slope = 1, linetype="dashed", color="grey") + facet_wrap(~genotype, scales="free") +
    labs(y="Estimated variance (95% CI)",x="Outcome variance conditional on genotype")
dev.off()

# absrsi modifier main effect with dummy method
results <- data.frame()
for (i in 1:n_sim){
    for (phi in seq(0, 6, 0.25)){
        theta <- phi * b
        x <- get_simulated_genotypes(q, n_obs)
        u <- rnorm(n_obs)
        y <- u*theta + x*u*theta + rnorm(n_obs)
        fit <- lm(y ~ x)
        d <- resid(fit)
        d2 <- abs(d)
        res <- data.frame(
            var0=var(y[x==0]),
            var1=var(y[x==1]),
            var2=var(y[x==2])
        )
        fit <- lm(d2 ~ as.factor(x))
        b0 <- fit %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(estimate)
        b1 <- fit %>% tidy %>% dplyr::filter(term == "as.factor(x)1") %>% dplyr::pull(estimate)
        b2 <- fit %>% tidy %>% dplyr::filter(term == "as.factor(x)2") %>% dplyr::pull(estimate)
        res$e0 <- b0^2/(2/pi)
        res$e1 <- b0^2/(2/pi) + (2*b0*b1+b1^2)/(2/pi)
        res$e2 <- b0^2/(2/pi) + (2*b0*b2+b2^2)/(2/pi)
        res$phi <- phi

        results <- rbind(results, res)
    }
}

e0 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e0) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==0)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
e1 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e1) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==1)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
e2 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(e2) %>% tidy %>% dplyr::mutate(genotype="var(Y|G==2)") %>% dplyr::select(estimate, conf.low, conf.high, genotype))
v0 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var0))
v1 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var1))
v2 <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var2))
e0 <- cbind(e0, v0 %>% dplyr::select(var))
e1 <- cbind(e1, v1 %>% dplyr::select(var))
e2 <- cbind(e2, v2 %>% dplyr::select(var))
e <- rbind(e0,e1,e2)

pdf("absrsi_dummy_model_modifier_main.pdf")
ggplot(e, aes(x=var, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() + geom_errorbar(width=.05) + theme_classic() + geom_abline(intercept = 0, slope = 1, linetype="dashed", color="grey") + facet_wrap(~genotype, scales="free") +
    labs(y="Estimated variance (95% CI)",x="Outcome variance conditional on genotype")
dev.off()