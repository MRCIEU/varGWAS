source("../../sim/funs.R")
set.seed(12345)

# QCTOOL on PATH

# sample size
n_obs <- 10000

# MAF
af <- 0.4

# covariates
data <- data.frame(
    S = paste0("S", seq(1, n_obs)),
    X = get_simulated_genotypes(af, n_obs),
    U = rnorm(n_obs),
    stringsAsFactors=F
)

# outcome
data$Y <- data$X * data$U + rnorm(n_obs)
data$Y <- data$Y / sd(data$Y)

# write out GEN file
write_gen("example.gen", "01", "SNPID_1", "RSID_1", "1", "A", "G", data$X)

# write phenotype & sample file
write.table(file = "example.csv", sep = ",", quote = F, row.names = F, data)

# convert to BGEN file & plink
system("qctool -g example.gen -og example.bgen")
system("bgenix -g example.bgen -clobber -index")