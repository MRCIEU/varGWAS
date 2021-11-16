library("data.table")
library("broom")
library("boot")
library('optparse')
source("funs.R")
#set.seed(12345)

option_list <- list(
  make_option(c("-p", "--phi"), type = "numeric", default = NULL, help = "Effect size of interaction relative to main effect"),
  make_option(c("-i", "--iteration"), type = "numeric", default = NULL, help = "Simulation iteration"),
  make_option(c("-n", "--n_iter"), type = "numeric", default = NULL, help = "Number of iterations within this sim")
);
opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

set.seed(opt$i + round(opt$p + 100))

# Requires OSCA and QCTOOL on PATH

# function to obtain regression weights
bs <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  result <- dummy_model(d$X, d$Y)
  return(result)
}

n_obs <- 1000
n_sim <- 200
af <- 0.4

# main effect size of X on Y detectable with 95% power
delta <- 0.169

# simulate GxE interaction effects and estimate power
results <- data.frame()
for (phi in seq(0, 6, 3)){
    theta <- delta * phi
    for (i in 1:n_sim) {
        message(paste0("phi:", phi, " i:", i))

        # simulate covariates
        data <- data.frame(
            S = paste0("S", seq(1, n_obs)),
            X = get_simulated_genotypes(af, n_obs),
            U = rnorm(n_obs),
            stringsAsFactors=F
        )

        # simulate outcome
        data$Y <- data$X * delta + data$U * delta + data$X * data$U * theta + rnorm(n_obs)
        data$Y <- scale(data$Y)

        # test for variance effect
        fit_boot <- boot(data=data, statistic=bs, R=50) %>% tidy
        fit_osca <- run_osca(data, T)
        fit_osca$BETA_x.osca_median <- fit_osca$BETA_x.osca_median / (2/pi)
        fit_osca$SE_x.osca_median <- fit_osca$SE_x.osca_median / (2/pi)
        #p <- lm(Y ~ X, data=data) %>% tidy %>% dplyr::filter(term == "X") %>% dplyr::pull(p.value)

        res <- data.frame(
            b0_dummy=fit_boot$statistic[1],
            s0_dummy=fit_boot$std.error[1],
            b1_dummy=fit_boot$statistic[2],
            s1_dummy=fit_boot$std.error[2],
            b2_dummy=fit_boot$statistic[3],
            s2_dummy=fit_boot$std.error[3],
            b1_osca=fit_osca$BETA_x.osca_median * 1,
            s1_osca=fit_osca$SE_x.osca_median * 1,
            b2_osca=fit_osca$BETA_x.osca_median * 2,
            s2_osca=fit_osca$SE_x.osca_median * 2
        )

        # add params
        #res <- data.frame(p)
        res$v0 <- var(data$Y[data$X==0])
        res$v1 <- var(data$Y[data$X==1])
        res$v2 <- var(data$Y[data$X==2])
        res$phi <- phi
        res$theta <- theta

        # store result
        results <- rbind(results, res)
    }
}

write.csv(results, file="sim12.csv")

get_est_ci <- function(results, estimand){
    results %>%
        dplyr::group_by(phi) %>%
        summarise_each_(
            funs_( 
                sprintf("t.test(%s) %>% tidy",estimand)
            ), 
            vars = vars_to_test)

        dplyr::summarize(t.test(paste0(estimand)) %>% tidy) %>%
        dplyr::select(phi, estimate, conf.low, conf.high) %>%
        dplyr::mutate(estimand=!!estimand)
}

b1_dummy_ci <- get_est_ci(results, "b1_dummy")



library("ggplot2")
ggplot(data=results, aes(x=lambda, y=est_power, ymin=est_power_low, ymax=est_power_high, group=phi, color=phi)) +
    geom_line() + 
    geom_point() + 
    geom_errorbar(width=.05) +
    theme_classic() + 
    xlab("Sample size inflation factor") + 
    ylab(paste0("Power (alpha=", 0.05, ")")) +
    labs(color = expression(phi)) +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks(n = 5)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey") +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey") +
    scale_color_viridis(direction = 1) +
    facet_grid(dist ~ method) +
    theme(
        strip.background = element_blank(),
        legend.title.align=0.5
    )