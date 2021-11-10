library("broom")
library("dplyr")
library("boot")
library('optparse')
source("funs.R")

option_list <- list(
  make_option(c("-b", "--beta"), type = "numeric", default = NULL, help = "Effect size of interaction"),
  make_option(c("-i", "--iteration"), type = "numeric", default = NULL, help = "Simulation iteration"),
  make_option(c("-n", "--n_iter"), type = "numeric", default = NULL, help = "Number of iterations within this sim")
);
opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

set.seed(opt$i + round(opt$b + 100))
n_obs <- 1000

# function to obtain regression weights
bs <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  result <- dummy_model(d$x, d$y)
  return(result)
}

df <- data.frame()
for (i in 1:opt$n){
    message(paste0("b:", opt$b, " i:", opt$i, " n:", opt$n, " j:", i))

    x <- get_simulated_genotypes(0.4, n_obs)
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

    df <- rbind(df, data.frame(b=opt$b, b1, lci1, uci1, b2, lci2, uci2, t1=var(y[x==1]) - var(y[x==0]), t2=var(y[x==2]) - var(y[x==0])))
}

write.table(df, file=paste0("results_i",opt$i,"_b",opt$b,".txt"))