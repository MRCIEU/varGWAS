library("data.table")
library("broom")
library("boot")
library('optparse')
source("funs.R")
set.seed(12345)

option_list <- list(
  make_option(c("-p", "--phi"), type = "numeric", default = NULL, help = "Effect size of interaction relative to main effect"),
  make_option(c("-i", "--int"), type = "numeric", default = NULL, help = "Linear=0, interaction=1")
);
opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

# Requires OSCA and QCTOOL on PATH

n_obs <- 10000
n_sim <- 1000
af <- 0.4

# main effect size of X on Y explaining 5% variance
delta <- 0.33
theta <- delta * opt$phi

# simulate GxE interaction effects and estimate power
results <- data.frame()
for (i in 1:n_sim) {
    message(paste0("phi:", opt$phi, " i:", i))

    # simulate covariates
    data <- data.frame(
        S = paste0("S", seq(1, n_obs)),
        X = get_simulated_genotypes(af, n_obs),
        U = rnorm(n_obs),
        stringsAsFactors=F
    )

    # simulate outcome
    if (opt$i == 1){
      data$Y <- data$X * delta + data$U * delta + data$X * data$U * theta + rnorm(n_obs)
    } else if (opt$i == 0){
      data$Y <- data$X * delta + data$U * delta + data$X * data$U * 0 + rnorm(n_obs, sd=sqrt(1 + data$X * theta))
    }
    data$Y <- scale(data$Y)

    # test for variance effect
    fit <- run_models(data)
    fit$BETA_x.osca_median <- fit$BETA_x.osca_median / (2/pi)
    fit$SE_x.osca_median <- fit$SE_x.osca_median / (2/pi)

    res <- data.frame(
        b1_dummy=fit$phi_x1,
        s1_dummy=fit$se_x1,
        b2_dummy=fit$phi_x2,
        s2_dummy=fit$se_x2,
        b1_osca=fit_osca$BETA_x.osca_median * 1,
        s1_osca=fit_osca$SE_x.osca_median * 1,
        b2_osca=fit_osca$BETA_x.osca_median * 2,
        s2_osca=fit_osca$SE_x.osca_median * 2
    )

    # add params
    res$v0 <- var(data$Y[data$X==0])
    res$v1 <- var(data$Y[data$X==1])
    res$v2 <- var(data$Y[data$X==2])
    res$phi <- opt$phi
    res$theta <- theta

    # store result
    results <- rbind(results, res)
}

write.csv(results, file=paste0("sim12_", opt$phi, "_", opt$i, ".csv"))