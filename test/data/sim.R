library("data.table")
library("broom")
library("pwr")
set.seed(12345)

n_obs <- 200
n_sim <- 20
alpha <- 0.05

#' Function to perform Breusch-Pagan test using t-test
#' @param x vector of genotype
#' @param y vector of response
bp_t <- function(x, y){
    fit1 <- lm(y ~ x)
    d <- resid(fit1)^2
    fit2 <- lm(d ~ x)
    fit2 <- tidy(fit2)
    return(data.frame(BETA.r=fit2$estimate[2], SE.r=fit2$std.error[2], P.r=fit2$p.value[2]))
}

#' Function to perform Breusch-Pagan test using f-test
#' @param x vector of genotype
#' @param y vector of response
bp_f <- function(x, y){
    fit1 <- lm(y ~ x)
    d <- resid(fit1)^2
    fit2 <- lm(d ~ x)
    fit0 <- lm(d ~ 1)
    f <- anova(fit0, fit2)
    fit2 <- tidy(fit2)
    f <- tidy(f)
    return(data.frame(BETA.r=fit2$estimate[2], SE.r=fit2$std.error[2], P.r=f$p.value[2]))
}

#' Function to simulate genotypes in HWE
#' @param q Recessive/alternative allele frequency
#' @param n_obs Number of observations to return
get_simulated_genotypes <- function(q, n_obs){
    p <- 1 - q
    x <- sample(c(0, 1, 2), n_obs, prob=c(p^2, 2 * p * q, q^2), replace=T)
    return(x)
}

# effect size to detect with assuming power, alpha and n_obs params
# TODO add AF
delta <- sqrt(pwr.f2.test(u = 1, v = n_obs - 1 - 1, sig.level = alpha, power = 0.8)$f2)

# simulate GxE interaction effects and estimate power
results <- data.frame()
for (phi in seq(1)){
    theta <- delta * phi
    beta <- delta - theta
    for (af in c(0.4)){
        for (lambda in c(1)){
            for (i in 1:n_sim){
                # simulate data
                x <- get_simulated_genotypes(af, n_obs * lambda)
                u <- rnorm(n_obs * lambda)
                y <- x*beta + x*u*theta + rnorm(n_obs * lambda)
                s <- paste0("S", seq(1, n_obs * lambda))
                A <- 0*0*1+1
                B <- 2*0*theta*1
                C <- theta*theta*1

                # write out GEN file
                fileConn<-file("genotypes.gen")
                writeLines(c(paste("01","SNPID_1", "RSID_1", "1", "A", "G", paste(sapply(x, function(g) if (g==0) { "0 0 0" } else if (g==1) {"0 1 0"} else if (g==2){"0 0 1"}), collapse=" "), collapse=" ")), fileConn)
                close(fileConn)
                write.csv(data.frame(s, x), file="genotypes.csv", quote=F, row.names=F)

                # convert to BGEN file
                system("qctool -g genotypes.gen -og genotypes.bgen 2> /dev/null")
                system("../../lib/bgen/build/apps/bgenix -g genotypes.bgen -clobber -index 2> /dev/null")
                
                # write phenotype file
                write.table(file="phenotypes.csv", sep=",", quote=F, row.names=F, data.frame(s, y))

                # run vGWAS using C++
                system("../../build/src/jlst_cpp_run -v phenotypes.csv -s , -o gwas.txt -b genotypes.bgen -p y -i s")

                # parse output
                res <- fread("gwas.txt", select=c("BETA", "SE", "P"), col.names=c("BETA.cpp", "SE.cpp", "P.cpp"))

                # run B-P using R
                res <- cbind(res, bp_t(x, y))

                # add params
                res$A<-A
                res$B<-B
                res$C<-C
                res$phi <- phi
                res$af <- af
                res$lambda <- lambda
                res$theta <- theta
                res$beta <- beta
                results <- rbind(results, res)
            }
        }
    }
}