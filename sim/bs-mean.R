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
#n_sim <- 200
#phi <- 2
#t1 <- 4
#t2 <- 16

# LAD-BF variance effects
model <- function(x, y, covar=NULL){
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
    # extract coef
    b0 <- fit2 %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull("estimate")
    b1 <- fit2 %>% tidy %>% dplyr::filter(term == "x1") %>% dplyr::pull("estimate")
    b2 <- fit2 %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull("estimate")
    # variance betas
    return(c(
        (2*b0*b1+b1^2)/(2/pi), # SNP=1
        (2*b0*b2+b2^2)/(2/pi) # SNP=2
    ))
}

# function to obtain regression weights
bs <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  result <- model(d$x, d$y)
  return(result)
}

df <- data.frame()
for (i in 1:opt$n){
    message(paste0("b:", opt$b, " i:", opt$i, " n:", opt$n, " j:", i))

    x <- rbinom(n_obs, 2, .5)
    u <- rnorm(n_obs)
    y <- x*u*opt$b + rnorm(n_obs)
    data <- data.frame(x,y)

    # bootstrapping with 1000 replications
    results <- boot(data=data, statistic=bs, R=1000)

    # get 95% confidence intervals
    ci1 <- boot.ci(results, type="bca", index=1)
    ci2 <- boot.ci(results, type="bca", index=2)

    # estimates
    b1 <- as.numeric(ci1$t0)
    lci1 <- ci1$bca[4]
    uci1 <- ci1$bca[5]
    b2 <- as.numeric(ci2$t0)
    lci2 <- ci2$bca[4]
    uci2 <- ci2$bca[5]

    df <- rbind(df, data.frame(b1, lci1, uci1, b2, lci2, uci2, t1=var(y[x==1]) - var(y[x==0]), t1=var(y[x==2]) - var(y[x==0])))
}

write.table(df, file=paste0("results_i",opt$i,"_b",opt$b,".txt"))

# check the coverage probability
#binom.test(sum(df$lci1 <= t1 & df$uci1 >= t1), n_sim)
#binom.test(sum(df$lci2 <= t2 & df$uci2 >= t2), n_sim)