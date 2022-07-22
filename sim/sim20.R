library("dplyr")
library("broom")
library("ggplot2")
library("varGWASR")
source("funs.R")
set.seed(123)

n_sim <- 200

results <- data.frame()
for (n_obs in c(1000, 5000, 25000)){
    for (n_mod in seq(1, 5)){
        p_lm <- c()
        p_bf <- c()
        for (i in 1:n_sim){
            # simulate phenotype
            x <- get_simulated_genotypes(0.4, n_obs)
            u <- matrix(rnorm(n_obs*n_mod), n_obs, n_mod)
            y <- x*.1275 + u[,1]*.135 + x*u[,1]*.13 + rnorm(n_obs)
            # test for interaction
            p_lm <- c(p_lm, lm(y ~ x*u) %>% tidy %>% dplyr::filter(grepl("x:u", term)) %>% dplyr::pull(p.value))
            p_bf <- c(p_bf, varGWASR::model(data.frame(x, y), "x", "y", covar1 = NULL, covar2 = NULL)$phi_p)
        }
        if (n_mod == 1){
            r_bf <- binom.test(sum(p_bf < 0.05), n_sim) %>% tidy %>% dplyr::mutate(n_obs=n_obs, n_mod=0, test="LAD-BF")
            results <- rbind(results, r_bf)
        }
        r_lm <- binom.test(sum(p_lm < 0.05), n_sim * n_mod) %>% tidy %>% dplyr::mutate(n_obs=n_obs, n_mod=n_mod, test="Linear regression")
        results <- rbind(results, r_lm)
    }
}

# plot
results$n_mod <- as.factor(results$n_mod)
pdf("sim20.pdf")
ggplot(results, aes(x=n_mod, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() +
    geom_errorbar(width=.05) +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks(n = 5)) +
    facet_grid(n_obs~test, scales="free_x", space="free_x") +
    xlab("Number of modifiers tested") +
    ylab("Power (95% CI)") +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey")
dev.off()