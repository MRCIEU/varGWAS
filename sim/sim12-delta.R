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

option_list <- list(
  make_option(c("-p", "--phi"), type = "numeric", default = NULL, help = "Effect size of interaction relative to main effect")
);
opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

get_se <- function(grad, vb){
    vG <- t(grad) %*% vb %*% grad
    sqrt(vG)
}

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
    b0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull("estimate")
    b1 <- fit2 %>% tidy %>% dplyr::filter(term == "x1") %>% dplyr::pull("estimate")
    b2 <- fit2 %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull("estimate")
    s1_dummy <- get_se(c(
        0,
        (2 * b0 + 2 * b1)/(2/pi),
        0
    ), vcovHC(fit2, type = "HC3"))
    s2_dummy <- get_se(c(
        0,
        0,
        (2 * b0 + 2 * b2)/(2/pi)
    ), vcovHC(fit2, type = "HC3"))
    return(data.frame(
      b1_dummy=(2*b0*b1+b1^2)/(2/pi),
      s1_dummy,
      b2_dummy=(2*b0*b2+b2^2)/(2/pi),
      s2_dummy
    ))
}

n_obs <- 100000
n_sim <- 2000
af <- 0.4

# main effect size of X on Y detectable with 95% power
delta <- 0.33

# simulate GxE interaction effects and estimate power
results <- data.frame()
#for (phi in seq(0, 6, .5)){
  theta <- delta * opt$phi

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
    data$Y <- data$X * delta + data$U * delta + data$X * data$U * theta + rnorm(n_obs)
    #data$Y <- scale(data$Y)

    # test for variance effect
    res <- dummy_model_delta(data$X, data$Y)

    # add params
    #res$v0 <- var(data$Y[data$X==0])
    res$v1 <- var(data$Y[data$X==1]) - var(data$Y[data$X==0])
    res$v2 <- var(data$Y[data$X==2]) - var(data$Y[data$X==0])
    res$phi <- opt$phi
    res$theta <- theta

    # store result
    results <- rbind(results, res)
  }

#}

write.csv(results, file=paste0("sim12-delta_", opt$phi, ".csv"))