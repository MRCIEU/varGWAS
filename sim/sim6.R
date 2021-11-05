library("broom")
library("cqrReg")
library("dplyr")
library("boot")
library('optparse')
library("car")
set.seed(23)

option_list <- list(
  make_option(c("-b", "--beta"), type = "numeric", default = NULL, help = "Effect size of interaction"),
  make_option(c("-i", "--iteration"), type = "numeric", default = NULL, help = "Simulation iteration"),
  make_option(c("-n", "--n_iter"), type = "numeric", default = NULL, help = "Number of iterations within this sim")
);
opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

n_obs <- 1000

get_residual <- function(x, y, covar=NULL){
    if (!is.null(covar)){
        X <- as.matrix(cbind(x, covar))
    } else {
        X <- as.matrix(data.frame(x))
    }
    # betas
    fit <- qrfit(X=X, y=y, tau=.5, method="mm")
    b <- rbind(fit$b, fit$beta)
    # predicted
    X <- cbind(rep(1, nrow(X)), X)
    fitted <- X %*% b
    # residual
    d <- y - fitted
    return(as.vector(d))
}

model <- function(data, indices){
    # subset vectors
    data <- data[indices,]
    x <- data$x
    y <- data$y
    # absresi
    d <- abs(get_residual(x, y))
    # dummy SNP
    x <- as.factor(x)
    # second-stage model
    fit2 <- lm(d ~ x)
    # extract coef
    b0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull("estimate")
    b1 <- fit2 %>% tidy %>% dplyr::filter(term == "x1") %>% dplyr::pull("estimate")
    b2 <- fit2 %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull("estimate")
    # difference between 0 vs 1 and 0 vs 2
    return(c(
        (2*b0*b1+b1^2)/(2/pi),
        (2*b0*b2+b2^2)/(2/pi)
    ))
}

results <- data.frame()
for (j in 1:opt$n){
    message(paste0("b:", opt$b, " i:", opt$i, " n:", opt$n, " j:", opt$j))
    # SNP
    x <- rbinom(n_obs, 2, .5)
    # modifier
    u <- rnorm(n_obs)
    # outcome
    y <- x*u*opt$b + rnorm(n_obs)
    # estimate variance effects and boostrap CI
    bs <- boot(data.frame(x, y), model, R=1000, stype="i")
    ci1 <- boot.ci(bs, type="bca", index=1)
    ci2 <- boot.ci(bs, type="bca", index=2)
    # extract parameters
    e1 <- bs %>% tidy %>% dplyr::pull("statistic") %>% dplyr::nth(1)
    se1 <- bs %>% tidy %>% dplyr::pull("std.error") %>% dplyr::nth(1)
    e2 <- bs %>% tidy %>% dplyr::pull("statistic") %>% dplyr::nth(2)
    se2 <- bs %>% tidy %>% dplyr::pull("std.error") %>% dplyr::nth(1)
    lci1 <- ci1$bca[4]
    uci1 <- ci1$bca[5]
    lci2 <- ci2$bca[4]
    uci2 <- ci2$bca[5]
    # store results
    results <- rbind(results, data.frame(
        v1=var(y[x==1]) - var(y[x==0]), # true SNP=0 vs SNP=1 variance
        v2=var(y[x==2]) - var(y[x==0]), # true SNP=0 vs SNP=2 variance
        e1, e2, se1, se2, # estimated SNP=0 vs SNP=1 & 2 variance
        lci1, uci1,
        lci2, uci2,
        b=opt$b,
        j,
        i=opt$i
    ))
}

write.table(results, file=paste0("results_i",opt$i,"_b",opt$b,".txt"))