library("broom")
library("ggplot2")
library("dplyr")
source("funs.R")

n_sim <- 200
n_obs <- 1000
q <- 0.25 # MAF
b <- 0.145 # main effect size to have 80% power

# interaction on change in var(Y)
results <- data.frame()
for (i in 1:n_sim){
    for (phi in seq(0, 12, 0.5)){
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

# linear effect of X on var(Y)
results <- data.frame()
for (i in 1:n_sim){
    for (phi in seq(0, 12, 0.5)){
        theta <- phi * b
        x <- get_simulated_genotypes(q, n_obs)
        u <- rnorm(n_obs)
        y <- x*u*theta + rnorm(n_obs)
        v1 <- var(y[x==1]) - var(y[x==0])
        v2 <- var(y[x==2]) - var(y[x==0])
        results <- rbind(results, data.frame(
            phi, v1, v2
        ))
    }
}
summary1 <- results %>%
    dplyr::group_by(phi) %>% 
    dplyr::summarise(t.test(v1) %>% tidy)
summary1$genotype <- "var(Y|G==1)"
summary2 <- results %>%
    dplyr::group_by(phi) %>% 
    dplyr::summarise(t.test(v2) %>% tidy)
summary2$genotype <- "var(Y|G==2)"
summary <- rbind(summary1,summary2)

pdf("Interaction_variance_gt_cond.pdf")
ggplot(summary, aes(x=phi, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() + geom_errorbar(width=.05) + theme_classic() + facet_grid(~genotype, scales="free_y")
    labs(y="Change in total variance (95% CI)",x="Interaction effect size relative to main effect",shape="Modifier main effect")
dev.off()