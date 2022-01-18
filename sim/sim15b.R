library("broom")
library("ggplot2")
library("dplyr")
source("funs.R")

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
    dplyr::summarise(t.test(results$v0) %>% tidy %>% dplyr::mutate(genotype="0"))
estimates1 <- results %>% 
    dplyr::group_by(model) %>% 
    dplyr::summarise(t.test(results$v1) %>% tidy %>% dplyr::mutate(genotype="1"))
estimates2 <- results %>% 
    dplyr::group_by(model) %>% 
    dplyr::summarise(t.test(results$v2) %>% tidy %>% dplyr::mutate(genotype="2"))
estimates <- rbind(
    estimates0,estimates1,estimates2
)
estimates$model <- factor(estimates$model, levels=c("XU", "X + XU", "U + XU"))

pdf("Genotype_variance.pdf")
ggplot(estimates, aes(x=genotype, y=estimate, ymin=conf.low, ymax=conf.high, group=model, shape=model)) +
    geom_point(position=position_dodge(width=0.3)) + geom_errorbar(width=.05) + theme_classic() + 
    labs(y="Variance (95% CI)",x="Genotype",shape="Model")
dev.off()

# define variance model
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
        fit <- lm(d2 ~ x2) %>% tidy %>% dplyr::filter(term == "x2")
        fit$var1 <- var(y[x==1]) - var(y[x==0])
        fit$var2 <- var(y[x==2]) - var(y[x==0])
        fit$phi <- phi
        results <- rbind(results, fit)
    }
}

estimate <- results %>% dplyr::group_by(phi) %>%
    dplyr::summarise(t.test(estimate) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::mutate(genotype="1"))
vestimate <- results %>% dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var1)) %>% dplyr::select(var)
vestimate2 <- results %>% dplyr::group_by(phi) %>%
    dplyr::summarise(var=mean(var2))
estimate <- cbind(estimate, vestimate)
estimate2 <- estimate
estimate2$estimate <- estimate2$estimate * 2
estimate2$conf.low <- estimate2$conf.low * 2
estimate2$conf.high <- estimate2$conf.high * 2
estimate2$var <- vestimate2$var
estimate2$genotype <- "2"
estimate <- rbind(estimate, estimate2)

pdf("Quadratic_model.pdf")
ggplot(estimate, aes(x=var, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point(position=position_dodge(width=0.3)) + geom_errorbar(width=.05) + theme_classic() + facet_wrap(~genotype, scales="free")
    labs(y="Variance (95% CI)",x="Genotype",shape="Model")
dev.off()