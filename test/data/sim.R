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
    beta <- delta - theta
    for (af in c(0.05, 0.1, 0.2, 0.4)){
        for (lambda in c(1, 10, 100, 1000, 10000)){
            for (i in 1:n_sim){
                # simulate data
                x <- get_simulated_genotypes(af, n_obs * lambda)
                u <- rnorm(n_obs * lambda)
                y <- x*beta + x*u*theta + rnorm(n_obs * lambda)
                s <- paste0("S", seq(1, n_obs * lambda))

                # write out GEN file
                fileConn<-file("genotypes.gen")
                writeLines(c(paste("01","rs123", "1", "A", "G", paste(sapply(x, function(g) if (g==0) { "0 0 0" } else if (g==1) {"0 1 0"} else if (g==2){"0 0 1"}), collapse=" "), collapse=" ")), fileConn)
                close(fileConn)

                # convert to BGEN file
                system("qctool -g genotypes.gen -og genotypes.bgen")
                system("../../lib/bgen/build/apps/bgenix -g genotypes.bgen -index")
                
                # write phenotype file
                write.table(file="phenotypes.csv", sep=",", quote=F, data.frame(s, y))

                # run vGWAS
                system("../../build/jlst_cpp -v phenotypes.csv -s , -o gwas.txt -p y -i s")

                # parse output
                res <- fread("gwas.txt")
                break
            }
        }
    }
}