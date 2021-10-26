library("data.table")
library("broom")
library("genpwr")
library('optparse')
library("jlst")
source("funs.R")
set.seed(12345)

# Requires OSCA and QCTOOL on PATH
n_obs <- 200
n_sim <- 200
af <- 0.4
phi <- 2
lambda <- 100

# main effect size of X on Y detectable with 80% power
delta <- as.numeric(genpwr.calc(calc = "es", model = "linear", ge.interaction = NULL, N = n_obs, k = NULL, MAF = af, Power = 0.8, Alpha = 0.05, sd_y = 1, True.Model = "Additive", Test.Model = "Additive")$ES_at_Alpha_0.05)

# simulate GxE interaction effects and estimate power
results <- data.frame()
theta <- delta * phi
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
    data$X * data$U * theta + + rnorm(n_obs * lambda)

    # run models
    res <- run_models(data)

    # add expected variance parameters
    res$EXP_x <- 2 * delta * theta
    res$EXP_xsq <- theta * theta

    # add params
    res$phi <- phi
    res$af <- af
    res$lambda <- lambda
    res$theta <- theta
    res$delta <- delta

    # store result
    results <- rbind(results, res)
}