library("data.table")
library("broom")
library("dplyr")
library("genpwr")
library("ggplot2")
source("funs.R")
set.seed(12345)

n_obs <- 200
n_sim <- 200
af <- 0.4

# simulate interaction effects and test using exhausitve pairwise analysis
results <- data.frame()
for (lambda in c(1, 10, 100, 1000)) {
  for (n_mod in c(1, 5, 10, 20)) {
    res <- data.frame()

    for (i in 1:n_sim) {
        # simulate modifiers that have no effect on outcome
        modifiers <- matrix(rnorm((n_mod-1)*n_obs*lambda), n_obs*lambda, n_mod-1)

        # simulate covariates
        data <- data.frame(
            S = paste0("S", seq(1, n_obs * lambda)),
            X = get_simulated_genotypes(af, n_obs * lambda),
            U = rnorm(n_obs * lambda),
            stringsAsFactors=F
        )

        # simulate outcome
        # interaction term to have 80% power when lambda = 1x
        data$Y <- data$X * data$U * 0.29 + rnorm(n_obs * lambda)

        # interaction test
        lm.fit <- lm(Y ~ X * U, data)
        lm_p <- apply(modifiers, 2, function(U) lm(data$Y ~ data$X * U) %>% tidy %>% dplyr::filter(term == "data$X:U") %>% dplyr::pull(p.value))

        # variance test
        phi_p <- rep(NA, n_mod)
        phi_p[1] <- varGWASR::model(data, "X", "Y")$phi_p

        res <- rbind(res,
            data.frame(
                lm_p=c(lm.fit %>% tidy %>% dplyr::filter(term == "X:U") %>% dplyr::pull(p.value), lm_p),
                phi_p,
                lambda,
                n_mod
            )
        )
    }
    
    # calculate power
    res$lm_h1 <- as.numeric(res$lm_p < 0.05)
    res$phi_h1 <- as.numeric(res$phi_p < 0.05)
    lm_pr <- tidy(binom.test(sum(res$lm_h1), nrow(res))) %>% dplyr::select(estimate, conf.low, conf.high)
    phi_pr <- tidy(binom.test(sum(res$phi_h1, na.rm=T), sum(!is.na(res$phi_p)))) %>% dplyr::select(estimate, conf.low, conf.high)
    lm_pr$n_mod <- n_mod
    lm_pr$lambda <- lambda
    lm_pr$model <- "Interaction test"
    phi_pr$n_mod <- n_mod
    phi_pr$lambda <- lambda
    phi_pr$model <- "Variance test"

    # store results
    results <- rbind(results, lm_pr)
    results <- rbind(results, phi_pr)
  }
}

# plot results
results$lambda <- as.factor(results$lambda)
ggplot(data=results, aes(x=lambda, y=estimate, ymin=conf.low, ymax=conf.high, group=model)) +
    geom_line() + 
    geom_point() + 
    geom_errorbar(width=.05) +
    theme_classic() + 
    xlab("Sample size inflation factor") + 
    ylab(paste0("Power (alpha=", 0.05, ")")) +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks(n = 5)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey") +
    facet_grid(n_mod ~ model) +
    theme(
        strip.background = element_blank(),
        legend.title.align=0.5
    )