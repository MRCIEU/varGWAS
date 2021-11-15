library("data.table")
library("broom")
library("boot")
source("funs.R")
set.seed(12345)

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

# main effect size of X on Y detectable with 80% power
delta <- 0.169

# simulate GxE interaction effects and estimate power
results <- data.frame()
for (phi in seq(6, 6, 0)){
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
        fit_boot <- boot(data=data, statistic=bs, R=75) %>% tidy
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

        # store result
        results <- rbind(results, res)
    }
}