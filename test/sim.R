library("data.table")
library("broom")
library("genpwr")
set.seed(12345)

n_obs <- 200
n_sim <- 1
af <- 0.4

#' Function to perform Breusch-Pagan test using t-test
#' @param x vector of genotype
#' @param y vector of response
bp_t <- function(x, y) {
  fit1 <- lm(y ~ x)
  d <- resid(fit1)^2
  fit2 <- lm(d ~ x)
  fit2 <- tidy(fit2)
  return(data.frame(BETA.r = fit2$estimate[2], SE.r = fit2$std.error[2], P.r = fit2$p.value[2]))
}

#' Function to perform Breusch-Pagan test using f-test
#' @param x vector of genotype
#' @param y vector of response
bp_f <- function(x, y) {
  fit1 <- lm(y ~ x)
  d <- resid(fit1)^2
  fit2 <- lm(d ~ x)
  fit0 <- lm(d ~ 1)
  f <- anova(fit0, fit2)
  fit2 <- tidy(fit2)
  f <- tidy(f)
  return(data.frame(BETA.r = fit2$estimate[2], SE.r = fit2$std.error[2], P.r = f$p.value[2]))
}

#' Function to simulate genotypes in HWE
#' @param q Recessive/alternative allele frequency
#' @param n_obs Number of observations to return
get_simulated_genotypes <- function(q, n_obs) {
  p <- 1 - q
  x <- sample(c(0, 1, 2), n_obs, prob = c(p^2, 2 * p * q, q^2), replace = T)
  return(x)
}

# main effect size of X on Y detectable with 80% power
delta <- as.numeric(genpwr.calc(calc = "es", model = "linear", ge.interaction = NULL, N = n_obs, k = NULL, MAF = af, Power = 0.8, Alpha = 0.05, sd_y = 1, True.Model = "Additive", Test.Model = "Additive")$ES_at_Alpha_0.05)

# simulate GxE interaction effects and estimate power
results <- data.frame()
for (phi in seq(0, 6, 0.5)) {
  theta <- delta * phi

  for (lambda in c(1, 10, 100, 1000)) {
    for (i in 1:n_sim) {
      # simulate data
      data <- data.frame(
        S = paste0("S", seq(1, n_obs * lambda)),
        X = get_simulated_genotypes(af, n_obs * lambda),
        U = rnorm(n_obs * lambda),
        sex = rbinom(n_obs * lambda, 1, .5),
        age = sample(30:70, n_obs * lambda, replace = T),
        PC = sapply(1:10, function(x) rnorm(n_obs * lambda))
      )
      data$Y <- data$X * delta + data$X * data$U * theta + rnorm(n_obs * lambda)

      # write out GEN file
      fileConn <- file("data/genotypes.gen")
      writeLines(c(paste("01", "SNPID_1", "RSID_1", "1", "A", "G", paste(sapply(data$X, function(g) if (g == 0) { "1 0 0" } else if (g == 1) { "0 1 0" } else if (g == 2) { "0 0 1" }), collapse = " "), collapse = " ")), fileConn)
      close(fileConn)

      # convert to BGEN file
      system("qctool -g data/genotypes.gen -og data/genotypes.bgen")
      system("bgenix -g data/genotypes.bgen -clobber -index")

      # write phenotype file
      write.table(file = "data/phenotypes.csv", sep = ",", quote = F, row.names = F, data)

      # run vGWAS using C++
      system("../build/bin/jlst_cpp -v data/phenotypes.csv -s , -o data/gwas.txt -b data/genotypes.bgen -p Y -i S -t 1")

      # parse output
      res <- fread("data/gwas.txt", select = c("BETA", "SE", "P"), col.names = c("BETA.cpp", "SE.cpp", "P.cpp"))

      # run B-P using R
      res <- cbind(res, bp_t(data$X, data$Y))

      # run LM
      fit <- tidy(lm(Y ~ X * U, data=data))
      res$BETA.x <- fit$estimate[2]
      res$BETA.u <- fit$estimate[3]
      res$BETA.xu <- fit$estimate[4]
      res$SE.x <- fit$std.error[2]
      res$SE.u <- fit$std.error[3]
      res$SE.xu <- fit$std.error[4]

      # add params
      res$phi <- phi
      res$af <- af
      res$lambda <- lambda
      res$theta <- theta
      res$delta <- delta
      res$var <- 2 * theta^2

      # store result
      results <- rbind(results, res)
    }
  }

}

# save results for plotting
write.csv(results, file = "data/results.csv")