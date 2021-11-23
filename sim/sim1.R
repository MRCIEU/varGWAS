library("data.table")
library("broom")
library("genpwr")
library('optparse')
library("jlst")
source("funs.R")
set.seed(12345)

# Requires OSCA and QCTOOL on PATH

option_list <- list(
  make_option(c("-d", "--dist"), type = "character", default = "n", help = "Outcome distribution", metavar = "character")
);
opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

n_obs <- 200
n_sim <- 200
af <- 0.4

# main effect size of X on Y detectable with 80% power
delta <- as.numeric(genpwr.calc(calc = "es", model = "linear", ge.interaction = NULL, N = n_obs, k = NULL, MAF = af, Power = 0.8, Alpha = 0.05, sd_y = 1, True.Model = "Additive", Test.Model = "Additive")$ES_at_Alpha_0.05)

# simulate GxE interaction effects and estimate power
results <- data.frame()
for (phi in seq(0, 6, 0.5)) {
  theta <- delta * phi

  for (lambda in c(1, 10, 100, 1000)) {
    for (i in 1:n_sim) {

      # simulate covariates
      data <- data.frame(
        S = paste0("S", seq(1, n_obs * lambda)),
        X = get_simulated_genotypes(af, n_obs * lambda),
        U = rnorm(n_obs * lambda),
        stringsAsFactors=F
      )

      # simulate outcome
      data$Y <- data$X * delta +
        data$U * delta +
        data$X * data$U * theta

      # add error term
      if (opt$dist == "n") {
        data$Y <- data$Y + rnorm(n_obs * lambda)
      } else if (opt$dist == "mn") {
        data$Y <- data$Y + c(rnorm(n_obs * lambda *.9), rnorm(n_obs * lambda * .1, mean=5))
      } else if (opt$dist == "t") {
        data$Y <- data$Y + rt(n_obs * lambda, 4)
      } else if (opt$dist == "l") {
        data$Y <- data$Y + rlnorm(n_obs * lambda)
      } else {
        stop(paste0("Distribution not implemented: ", opt$dist))
      }

      # run models
      res <- run_models(data)

      # add params
      res$phi <- phi
      res$af <- af
      res$lambda <- lambda
      res$theta <- theta
      res$delta <- delta

      # store result
      results <- rbind(results, res)
    }
  }

}

# print warnings
warnings()

# save results for plotting
write.csv(results, file = paste0("data/power_", opt$dist, ".csv"))