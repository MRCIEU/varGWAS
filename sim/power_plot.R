library("ggplot2")
library("dplyr")
library("tidyr")
library("broom")
library('optparse')
library("data.table")
library("viridis")
set.seed(123)

#' Function to estimate the power of an MC experiment
#' @param results Dataframe containing: P value and analysis group(s)
#' @param field Name of P value to summarise
#' @param n_sim Number of simulations performed
#' @param grp_name A vector of fields for use in grouping the analysis
#' @param alpha P threshold
calc_power <- function(results, field, n_sim, grp_name, alpha=0.05){
    # threshold P value
    results <- results %>%
        mutate(pos = as.numeric(get(field) < alpha))
    
    # count positives and estimate power with 95% CI
    h1 <- results %>%
        select(pos, all_of(grp_name)) %>%
        drop_na() %>%
        group_by_at(vars(all_of(grp_name))) %>%
        summarise(h1 = sum(pos)) %>%
        rowwise() %>%
        mutate(est_power = as.numeric(tidy(binom.test(h1, n_sim))$estimate)) %>% 
        mutate(est_power_low = as.numeric(tidy(binom.test(h1, n_sim))$conf.low)) %>%
        mutate(est_power_high = as.numeric(tidy(binom.test(h1, n_sim))$conf.high))

    return(h1 %>% ungroup())
}

# load data
n <- fread(paste0("data/power_n.csv"))
l <- fread(paste0("data/power_l.csv"))
t <- fread(paste0("data/power_t.csv"))

# process data
cpp_bp_n <- calc_power(n, "P.cpp_bp", 200, c("phi", "lambda"))
cpp_bp_n$dist <- "Normal"
cpp_bp_n$method <- "Breusch-Pagan"
osca_n <- calc_power(n, "P.osca", 200, c("phi", "lambda"))
osca_n$dist <- "Normal"
osca_n$method <- "OSCA"
cpp_bp_l <- calc_power(l, "P.cpp_bp", 200, c("phi", "lambda"))
cpp_bp_l$dist <- "Lognormal"
cpp_bp_l$method <- "Breusch-Pagan"
osca_l <- calc_power(l, "P.osca", 200, c("phi", "lambda"))
osca_l$dist <- "Lognormal"
osca_l$method <- "OSCA"
cpp_bp_t <- calc_power(t, "P.cpp_bp", 200, c("phi", "lambda"))
cpp_bp_t$dist <- "T-dist"
cpp_bp_t$method <- "Breusch-Pagan"
osca_t <- calc_power(t, "P.osca", 200, c("phi", "lambda"))
osca_t$dist <- "T-dist"
osca_t$method <- "OSCA"
results <- rbind(cpp_bp_n, cpp_bp_l, cpp_bp_t, osca_n, osca_l, osca_t)
results$method <- factor(results$method, levels = c("Breusch-Pagan", "OSCA"))
results$dist <- factor(results$dist, levels = c("Normal", "Lognormal", "T-dist"))
results$lambda <- factor(results$lambda, levels = c(1,10,100,1000))

# create plot
pdf("power.pdf")
ggplot(data=results, aes(x=lambda, y=est_power, ymin=est_power_low, ymax=est_power_high, group=phi, color=phi)) +
    geom_line() + 
    geom_point() + 
    geom_errorbar(width=.05) +
    theme_classic() + 
    xlab("Sample size inflation factor") + 
    ylab(paste0("Power (alpha=", 0.05, ")")) +
    labs(color = expression(phi)) +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks(n = 10)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey") +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey") +
    scale_color_viridis(direction = 1) +
    facet_grid(dist ~ method) +
    theme(
        strip.background = element_blank(),
        legend.title.align=0.5
    )
dev.off()