library("data.table")
library("broom")
library("genpwr")
library("jlst")
source("funs.R")
set.seed(12345)

n_obs <- 200
n_sim <- 200
af <- 0.4

# main effect size of X on Y detectable with 80% power
delta_x <- as.numeric(genpwr.calc(calc = "es", model = "linear", ge.interaction = NULL, N = n_obs, k = NULL, MAF = af, Power = 0.8, Alpha = 0.05, sd_y = 1, True.Model = "Additive", Test.Model = "Additive")$ES_at_Alpha_0.05)

# main effect size of confounder detectable with 80% power
delta_c1 <- get_binary_delta(n_obs, .5, 1, .05, .8)
delta_c2 <- get_cont_delta(n_obs, .05, .8)

# simulate GxE interaction effects and estimate power
results <- data.frame()
for (phi in seq(0, 6, 0.5)) {
  theta <- delta_x * phi

  for (lambda in c(1, 10, 100, 1000)) {
    for (i in 1:n_sim) {

      # simulate confounders
      C1 <- rbinom(n_obs * lambda, 1, .5)
      C2 <- sample(30:70, n_obs * lambda, replace = T)

      # simulate covariates
      data <- data.frame(
        S = paste0("S", seq(1, n_obs * lambda)),
        X = get_simulated_genotypes(af, n_obs * lambda),
        U = C1 * delta_c1 + C2 * delta_c2 + rnorm(n_obs * lambda),
        C1, C2
      )

      # simulate outcome
      data$Y <-
        data$C1 * delta_c1 +
        data$C2 * delta_c2 +
        data$X * delta_x +
        data$U * delta_x +
        data$X * data$U * theta +
        rnorm(n_obs * lambda)

      # write out GEN file
      write_gen("data/genotypes.gen", "01", "SNPID_1", "RSID_1", "1", "A", "G", data$X)

      # write phenotype & sample file
      write.table(file = "data/phenotypes.csv", sep = ",", quote = F, row.names = F, data)

      # convert to BGEN file & plink
      system("qctool -g data/genotypes.gen -og data/genotypes.bgen")
      system("bgenix -g data/genotypes.bgen -clobber -index")

      # run vGWAS
      system("jlst_cpp -v data/phenotypes.csv -s , -o data/gwas.txt -b data/genotypes.bgen -p Y -i S -t 1 -c C1,C2")

      # B-P using R implementation
      bp <- vartest(data$Y, x=data$X, covar=data.frame(data$C1, data$C2), type=1, x.sq=T, covar.var=T)
      res_r <- data.frame(BETA_x.r = bp$coef[2,1], SE_x.r = bp$coef[2,2], BETA_xsq.r = bp$coef[3,1], SE_xsq.r = bp$coef[3,2], P.r = as.numeric(bp$test[3]))

      # B-P using C++
      res_cpp <- fread("data/gwas.txt", select = c("beta", "se", "p", "phi_x", "se_x", "phi_xsq", "se_xsq", "phi_p"), col.names = c("BETA_mu.cpp", "SE_mu.cpp", "P_mu.cpp", "BETA_x.cpp", "SE_x.cpp", "BETA_xsq.cpp", "SE_xsq.cpp", "P.cpp"))

      # combine results
      res <- cbind(res_r, res_cpp)

      # add LM
      fit <- tidy(lm(Y ~ X * U + C1 + C2, data=data))
      res$BETA.x <- fit$estimate[2]
      res$BETA.u <- fit$estimate[3]
      res$BETA.xu <- fit$estimate[4]
      res$BETA.c1 <- fit$estimate[5]
      res$BETA.c2 <- fit$estimate[6]
      res$SE.x <- fit$std.error[2]
      res$SE.u <- fit$std.error[3]
      res$SE.xu <- fit$std.error[4]
      res$SE.c1 <- fit$std.error[5]
      res$SE.c2 <- fit$std.error[6]
      fit <- tidy(lm(Y ~ X + C1 + C2, data=data))
      res$BETA_mu.r <- fit$estimate[2]
      res$SE_mu.r <- fit$std.error[2]
      res$P_mu.r <- fit$p.value[2]

      # add expected variance parameters
      res$EXP_x <- 2*delta_x*theta
      res$EXP_xsq <- theta*theta

      # add params
      res$phi <- phi
      res$af <- af
      res$lambda <- lambda
      res$theta <- theta
      res$delta_x <- delta_x

      # store result
      results <- rbind(results, res)
    }
  }

}

# save results for plotting
write.csv(results, file = paste0("data/power_", opt$dist, ".csv"))