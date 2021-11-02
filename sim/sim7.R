library("data.table")
library("jlst")
library("broom")
library("optparse")

option_list <- list(
  make_option(c("-t", "--trait"), type = "character", default = "n", help = "Trait", metavar = "character")
);
opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

# read in emperical distribution
d <- fread(paste0("data/", opt$t, ".txt"))
n_sim <- 1000
n_obs <- 100000
af <- 0.05

# no effect
p <- rep(NA, n_sim)
for (i in 1:n_sim){
    s <- d[sample(1:nrow(d), replace=T, size=n_obs)]
    x <- rbinom(n_obs, 2, af)
    p[i] <- vartest(s[[opt$t]], x, covar=s %>% dplyr::select("age_at_recruitment.21022.0.0", "sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), covar.var = T, type = 2, x.sq = T)$test$P
}
binom.test(sum(p<0.05), n_sim) %>% tidy