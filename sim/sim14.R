library("data.table")
library("broom")
library("dplyr")
library("genpwr")
source("funs.R")
set.seed(12345)

n_obs <- 200
n_sim <- 200
af <- 0.4

# main effect size of X on Y detectable with 80% power
delta <- as.numeric(genpwr.calc(calc = "es", model = "linear", ge.interaction = NULL, N = n_obs, k = NULL, MAF = af, Power = 0.8, Alpha = 0.05, sd_y = 1, True.Model = "Additive", Test.Model = "Additive")$ES_at_Alpha_0.05)

# size of interaction relative to main effect
phi <- 0.5

# interaction effect size
theta <- delta * phi

# number of modifiers to evaluate
n_mods <- c(1, 5, 10, 20)

# simulate interaction effects and test using exhausitve pairwise analysis
results <- data.frame()
for (lambda in c(5, 50, 500)) {
  for (n in n_mods) {
    res <- data.frame()

    for (i in 1:n_sim) {
        # simulate modifiers that have no effect on outcome
        modifiers <- matrix(rnorm((n-1)*n_obs*lambda), n_obs*lambda, n-1)

        # simulate covariates
        data <- data.frame(
            S = paste0("S", seq(1, n_obs * lambda)),
            X = get_simulated_genotypes(af, n_obs * lambda),
            U = rnorm(n_obs * lambda),
            stringsAsFactors=F
        )

        # simulate outcome
        data$Y <- data$X * delta + data$U * delta + data$X * data$U * theta + rnorm(n_obs * lambda)

        # interaction test
        lm.fit <- lm(Y ~ X * U, data)
        lm_p <- apply(modifiers, 2, function(U) lm(data$Y ~ data$X * U) %>% tidy %>% dplyr::filter(term == "data$X:U") %>% dplyr::pull(p.value))

        # variance test
        phi_p <- rep(NA, n)
        phi_p[1] <- varGWASR::model(data, "X", "Y")$phi_p

        res <- rbind(res,
            data.frame(
                lm_p=c(lm.fit %>% tidy %>% dplyr::filter(term == "X:U") %>% dplyr::pull(p.value), lm_p),
                phi_p,
                lambda,
                n
            )
        )
    }
    
    # calculate power
    results$lm_h1 <- as.numeric(results$lm_p < 0.05)
    results$phi_h1 <- as.numeric(results$phi_p < 0.05)
    lm_pr <- tidy(binom.test(sum(results$lm_h1), nrow(results))) %>% dplyr::select(estimate, conf.low, conf.high)
    phi_pr <- tidy(binom.test(sum(results$phi_h1, na.rm=T), sum(!is.na(results$phi_p)))) %>% dplyr::select(estimate, conf.low, conf.high)
    lm_pr$n_mod <- n
    lm_pr$model <- "LM"
    phi_pr$n_mod <- n
    phi_pr$model <- "LAD-BF"

    # store results
    results <- rbind(results, lm_pr)
    results <- rbind(results, phi_pr)
  }
}