library("dplyr")
library("broom")
library("tidyr")
library("ggpubr")
library("lmtest")
library("data.table")
source("funs.R")
set.seed(123)

# Requires OSCA and QCTOOL on PATH

# TODO mult SNP analysis
n_sim <- 200
n_obs <- 300000

d <- data.frame()
for (i in 1:n_sim){
    # simulate covariates
    data <- data.frame(
        S = paste0("S", seq(1, n_obs)),
        X = get_simulated_genotypes(.4, n_obs),
        stringsAsFactors=F
    )

    # simulate outcome
    data$Y <- data$X + rlnorm(n_obs)

    # run models
    res <- run_models(data)

    # store result
    d <- rbind(d, res)
}

# save data
write.table(d, "data/sim4.txt")

# mean and 95% CI of elapsed time
results <- data.frame()
results <- rbind(results, tidy(t.test(d$elapsed.cpp_bp)))
results <- rbind(results, tidy(t.test(d$elapsed.osca_median)))
results <- rbind(results, tidy(t.test(d$elapsed.osca_mean)))
results <- rbind(results, tidy(t.test(d$elapsed.cpp_bf)))
results$method <- factor(c("Breusch-Pagan", "Brown-Forsythe", "Levene", "Brown-Forsythe (LAD)"),
    levels=c("Levene", "Breusch-Pagan", "Brown-Forsythe", "Brown-Forsythe (LAD)"))
results$location <- factor(c("Mean", "Median", "Mean", "Median"), levels=c("Mean", "Median"))

# barchart
pdf("data/sim4a.pdf")
ggplot(data=results, aes(x=method, y=estimate, ymin=conf.low, ymax=conf.high, group=location, color=location)) +
    geom_point() + 
    geom_errorbar(width=.05) +
    theme_classic() + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    labs(color="Location") +
    xlab("Method") +
    ylab("Mean runtime (seconds, 95% CI)")
dev.off()