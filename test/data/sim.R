library("data.table")
library("pwr")
set.seed(12345)

n_obs <- 200
n_sim <- 200
alpha <- 0.05

#' Function to simulate genotypes in HWE
#' @param q Recessive/alternative allele frequency
#' @param n_obs Number of observations to return
get_simulated_genotypes <- function(q, n_obs){
    p <- 1 - q
    x <- sample(c(0, 1, 2), n_obs, prob=c(p^2, 2 * p * q, q^2), replace=T)
    return(x)
}

# effect size to detect with assuming power, alpha and n_obs params
delta <- sqrt(pwr.f2.test(u = 1, v = n_obs - 1 - 1, sig.level = alpha, power = 0.8)$f2)

# simulate GxE interaction effects and estimate power
results <- data.frame()
for (phi in seq(0, 6, 0.5)){
    theta <- delta * phi
    b1 <- delta - theta
    for (af in c(0.05, 0.1, 0.2, 0.4)){
        for (lambda in c(1, 10, 100, 1000, 10000)){
            for (i in 1:n_sim){
                # simulate data
                x <- get_simulated_genotypes(af, n_obs * lambda)
                u <- rnorm(n_obs * lambda)
                y <- x*b1 + x*u*theta + rnorm(n_obs * lambda)
                s <- paste0("S", seq(1, n_obs * lambda))

                # write out GEN file
                g=sapply(x, function(g) if (g==0) { "0 0 0" } else if (g==1) {"0 1 0"} else if (g==2){"0 0 1"})
                df=data.frame(chr=rep("01", n_obs * lambda), rsid=paste0("rs", seq(1, n_obs * lambda)), g)

                # write BGEN file
                system(paste0("qctool -g example_#.gen -og example.bgen"))
                
                # write phenotype file
                write.table(file=paste0(), data.frame(s, y))
            }
        }
    }
}

# run vGWAS
system(paste0("../../build/jlst_cpp"))