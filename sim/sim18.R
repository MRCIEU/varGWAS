library("broom")
library("dplyr")
library("tidyr")
set.seed(1234)

n_sim <- 200
n_obs <- 200
alpha <- 0.05
power <- 0.8
phi <- 0.5
lambda <- c(100, 500, 3500)
n_mod <- c(1, 5, 10, 20)
delta <- get_cont_delta(n_obs, alpha, power)

# plan
#1) one variable has a true interaction effect with exposure on outcome
#2) vary the number of extra variables tested, i.e. variables that don’t have an effect
#3) calculate the power of that approach per N (1,5,10,20)?
#4) compare to approach that says test for difference in variances, then only do the LM tests if that test is significant

#Want to compare not just power but type I error too.

#Then also see what this table looks like if the sample size inflation factor is 100 and 500 (because then everything should be underpowered)

# power func
power_func <- function(l, n){
    results <- data.frame()

    for (i in 1:n_sim){
        # simulate cont phenotype with x cont interaction
        pheno <- sim_phenotype(n_obs * l, delta, phi, 0, 0, bin_exp=FALSE)

        # simulate modifiers that have no effect on outcome
        modifiers <- matrix(rnorm((n-1)*n_obs*l), n_obs*l, n-1)

        # interaction test
        lm.fit <- lm(y ~ x * u, pheno)
        lm_p <- apply(modifiers, 2, function(u) tidy(lm(pheno$y ~ pheno$x * u))$p.value[4])

        # variance test
        bp_p <- rep(NA, n)
        bp_p[1] <- bp(pheno$x, pheno$y)$test$p.value[2]

        results <- rbind(results,
            data.frame(
                lm_p=c(tidy(lm.fit)$p.value[4], lm_p),
                bp_p,
                l,
                n
            )
        )
    }

    results$lm_h1 <- as.numeric(results$lm_p < alpha)
    results$bp_h1 <- as.numeric(results$bp_p < alpha)
    lm_pr <- tidy(binom.test(sum(results$lm_h1), nrow(results)))
    bp_pr <- tidy(binom.test(sum(results$bp_h1, na.rm=T), sum(!is.na(results$bp_p))))

    return(data.frame(l, n, lm_pwr=lm_pr$estimate, bp_pwr=bp_pr$estimate))
}

# simulate interaction effects and test using exhausitve pairwise analysis
results <- data.frame()
for (l in lambda){
    for (n in c(1, 5, 10, 20)){
        results <- rbind(results, power_func(l, n))
    }
}

# write results to table
write.table(results, file="power_lm_vs_bp.txt", sep="\t", quote=F, row.names=F)