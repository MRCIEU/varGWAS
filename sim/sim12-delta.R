library("data.table")
library("broom")
library("boot")
library('optparse')
library("lmtest")
library("sandwich")
library('car')
library('cqrReg')
source("funs.R")
set.seed(12345)

get_se <- function(grad, vb){
    vG <- t(grad) %*% vb %*% grad
    sqrt(vG)
}

# LAD-BF variance effects
dummy_model_delta <- function(x, y, covar=NULL){
    if (!is.null(covar)){
        X <- as.matrix(cbind(x, covar))
    } else {
        X <- as.matrix(data.frame(x))
    }
    # first-stage fit
    fit <- qrfit(X=X, y=y, tau=.5, method="mm")
    b <- rbind(fit$b, fit$beta)
    # predicted
    X <- cbind(rep(1, nrow(X)), X)
    fitted <- X %*% b
    # residual
    d <- y - fitted
    # abs residual
    d <- abs(as.vector(d))
    # dummy SNP
    x <- as.factor(x)
    # second-stage model
    fit2 <- lm(d ~ x)
    # delta
    b0 <- deltaMethod(fit2, "b0^2/(2/pi)", parameterNames=c("b0", "b1", "b2"), vcov=vcovHC(fit2, type = "HC0"))
    b1 <- deltaMethod(fit2, "b0^2/(2/pi) + (2*b0*b1+b1^2)/(2/pi)", parameterNames=c("b0", "b1", "b2"), vcov=vcovHC(fit2, type = "HC0"))
    b2 <- deltaMethod(fit2, "b0^2/(2/pi) + (2*b0*b2+b2^2)/(2/pi)", parameterNames=c("b0", "b1", "b2"), vcov=vcovHC(fit2, type = "HC0"))
    # variance betas
    return(data.frame(
      b0_dummy=b0$Estimate,
      s0_dummy=b0$SE,
      b1_dummy=b1$Estimate,
      s1_dummy=b1$SE,
      b2_dummy=b2$Estimate,
      s2_dummy=b2$SE
    ))
}
dummy_model_delta_manual <- function(x, y, covar=NULL){
    if (!is.null(covar)){
        X <- as.matrix(cbind(x, covar))
    } else {
        X <- as.matrix(data.frame(x))
    }
    # first-stage fit
    fit <- qrfit(X=X, y=y, tau=.5, method="mm")
    b <- rbind(fit$b, fit$beta)
    # predicted
    X <- cbind(rep(1, nrow(X)), X)
    fitted <- X %*% b
    # residual
    d <- y - fitted
    # abs residual
    d <- abs(as.vector(d))
    # dummy SNP
    x <- as.factor(x)
    # second-stage model
    fit2 <- lm(d ~ x)
    # delta
    b0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull("estimate")
    b1 <- fit2 %>% tidy %>% dplyr::filter(term == "x1") %>% dplyr::pull("estimate")
    b2 <- fit2 %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull("estimate")
    b0_dm <- deltaMethod(fit2, "b0^2/(2/pi)", parameterNames=c("b0", "b1", "b2"), vcov=vcovHC(fit2, type = "HC0"))
    b1_dm <- deltaMethod(fit2, "b0^2/(2/pi) + (2*b0*b1+b1^2)/(2/pi)", parameterNames=c("b0", "b1", "b2"), vcov=vcovHC(fit2, type = "HC0"))
    b2_dm <- deltaMethod(fit2, "b0^2/(2/pi) + (2*b0*b2+b2^2)/(2/pi)", parameterNames=c("b0", "b1", "b2"), vcov=vcovHC(fit2, type = "HC0"))
    s0_dummy <- get_se(c(
        2 * b0/(2/pi),
        0,
        0
    ), vcovHC(fit2, type = "HC0"))
    s1_dummy <- get_se(c(
        2 * b0/(2/pi),
        (2 * b0 + 2 * b1)/(2/pi),
        0
    ), vcovHC(fit2, type = "HC0"))
    s2_dummy <- get_se(c(
        2 * b0/(2/pi),
        0,
        (2 * b0 + 2 * b2)/(2/pi)
    ), vcovHC(fit2, type = "HC0"))
    # variance betas
    return(data.frame(
      b0_dummy=b0^2/(2/pi),
      s0_dummy,
      b1_dummy=b0^2/(2/pi) + (2*b0*b1+b1^2)/(2/pi),
      s1_dummy,
      b2_dummy=b0^2/(2/pi) + (2*b0*b2+b2^2)/(2/pi),
      s2_dummy,
      b0_dm=b0_dm$Estimate,
      s0_dm=b0_dm$SE,
      b1_dm=b1_dm$Estimate,
      s1_dm=b1_dm$SE,
      b2_dm=b2_dm$Estimate,
      s2_dm=b2_dm$SE
    ))
}

n_obs <- 1000
n_sim <- 200
af <- 0.4

# main effect size of X on Y detectable with 95% power
delta <- 0.169

# simulate GxE interaction effects and estimate power
results <- data.frame()
for (phi in seq(0, 6, .5)){
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
    res <- dummy_model_delta_manual(data$X, data$Y)

    # add params
    res$v0 <- var(data$Y[data$X==0])
    res$v1 <- var(data$Y[data$X==1])
    res$v2 <- var(data$Y[data$X==2])
    res$phi <- phi
    res$theta <- theta

    # store result
    results <- rbind(results, res)
  }

}

# mean variance effect
v0_mean <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarize(v0_mean=mean(v0))
v1_mean <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarize(v1_mean=mean(v1))
v2_mean <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarize(v2_mean=mean(v2))
results <- merge(results, v0_mean, "phi")
results <- merge(results, v1_mean, "phi")
results <- merge(results, v2_mean, "phi")

# CI
results$b0_dummy_lci <- results$b0_dummy - (results$s0_dummy * 1.96)
results$b0_dummy_uci <- results$b0_dummy + (results$s0_dummy * 1.96)
results$b1_dummy_lci <- results$b1_dummy - (results$s1_dummy * 1.96)
results$b1_dummy_uci <- results$b1_dummy + (results$s1_dummy * 1.96)
results$b2_dummy_lci <- results$b2_dummy - (results$s2_dummy * 1.96)
results$b2_dummy_uci <- results$b2_dummy + (results$s2_dummy * 1.96)

# coverage
b0_dummy <- results %>% 
    dplyr::group_by(phi) %>%
    dplyr::summarize(binom.test(sum(b0_dummy_lci <= v0_mean & b0_dummy_uci >= v0_mean), n()) %>% tidy) %>%
    dplyr::select(phi, estimate, conf.low, conf.high) %>%
    dplyr::mutate(genotype="SNP=0", method="Dummy-delta")

b1_dummy <- results %>% 
    dplyr::group_by(phi) %>%
    dplyr::summarize(binom.test(sum(b1_dummy_lci <= v1_mean & b1_dummy_uci >= v1_mean), n()) %>% tidy) %>%
    dplyr::select(phi, estimate, conf.low, conf.high) %>%
    dplyr::mutate(genotype="SNP=1", method="Dummy-delta")

b2_dummy <- results %>% 
    dplyr::group_by(phi) %>%
    dplyr::summarize(binom.test(sum(b2_dummy_lci <= v2_mean & b2_dummy_uci >= v2_mean), n()) %>% tidy) %>%
    dplyr::select(phi, estimate, conf.low, conf.high) %>%
    dplyr::mutate(genotype="SNP=2", method="Dummy-delta")
  
coverage <- rbind(b0_dummy, b1_dummy, b2_dummy)

library("ggplot2")
ggplot(data=coverage, aes(x=phi, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() + 
    geom_errorbar(width=.05) +
    theme_classic() +
    xlab("Phi") +
    ylab("Coverage of 95% CI (95% CI)") +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey") +
    facet_grid(method ~ genotype, scales="free_x") +
    labs(shape="Genotype") +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks(n = 5)) +
    theme(
        strip.background = element_blank(),
        legend.title.align=0.5
    )