library("data.table")
library("dplyr")
library("ggplot2")
library("broom")
library("scales")
set.seed(1234)

# load simulation data
d <- fread("data/sim2.csv")

# mean and 95% CI of elapsed time
results <- data.frame()
results <- rbind(results, tidy(t.test(d$elapsed.cpp_bp)))
results <- rbind(results, tidy(t.test(d$elapsed.osca_median)))
results <- rbind(results, tidy(t.test(d$elapsed.osca_mean)))
results <- rbind(results, tidy(t.test(d$elapsed.cpp_bf)))
results$method <- factor(c("Levene(OLS)", "Levene(median)", "Levene(mean)", "Levene(LAD)"),
    levels=c("Levene(mean)", "Levene(OLS)", "Levene(median)", "Levene(LAD)"))
results$residual <- factor(c("Mean", "Median", "Mean", "Median"), levels=c("Mean", "Median"))

# barchart
ggplot(data=results, aes(x=method, y=estimate, ymin=conf.low, ymax=conf.high, group=residual, color=residual)) +
    geom_point() + 
    geom_errorbar(width=.05) +
    theme_classic() + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    labs(color="Residual") +
    xlab("Method") +
    ylab("Mean runtime (seconds, 95% CI)")