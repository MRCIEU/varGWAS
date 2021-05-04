library("data.table")
library("broom")
library("genpwr")
library('optparse')
set.seed(12345)

option_list <- list(
  make_option(c("-d", "--dist"), type="character", default=NULL, help="Outcome distribution", metavar="character")
);
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

n_obs <- 200
n_sim <- 200
af <- 0.4

#' Function to perform Breusch-Pagan test using f-test
#' @param x vector of genotype
#' @param y vector of response
bp <- function(x, y) {
  xsq <- x^2
  fit1 <- lm(y ~ x)
  d <- resid(fit1)^2
  fit2 <- lm(d ~ x + xsq)
  fit0 <- lm(d ~ 1)
  f <- anova(fit0, fit2)
  fit2 <- tidy(fit2)
  f <- tidy(f)
  return(data.frame(BETA_x.r = fit2$estimate[2], SE_x.r = fit2$std.error[2], BETA_xsq.r = fit2$estimate[3], SE_xsq.r = fit2$std.error[3], P.r = f$p.value[2]))
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

      if (opt$dist == "n"){
        data$Y <- data$X * delta + data$U * delta + data$X * data$U * theta + rnorm(n_obs * lambda)
      } else if (opt$dist == "t"){
        data$Y <- data$X * delta + data$U * delta + data$X * data$U * theta + rt(n_obs * lambda, 4)
      } else if (opt$dist == "l"){
        data$Y <- data$X * delta + data$U * delta + data$X * data$U * theta + rlnorm(n_obs * lambda)
      } else {
        stop(paste0("Distribution not implemented: ", opt$dist))
      }

      # write out GEN file
      fileConn <- file("data/genotypes.gen")
      writeLines(c(paste("01", "SNPID_1", "RSID_1", "1", "A", "G", paste(sapply(data$X, function(g) if (g == 0) { "1 0 0" } else if (g == 1) { "0 1 0" } else if (g == 2) { "0 0 1" }), collapse = " "), collapse = " ")), fileConn)
      close(fileConn)

      # write phenotype & sample file
      write.table(file = "data/phenotypes.csv", sep = ",", quote = F, row.names = F, data)
      write.table(file = "data/phenotypes.txt", sep = "\t", quote = F, row.names = F, col.names = F, data[,c("S", "S", "Y")])
      fileConn <- file("data/samples.txt")
      writeLines(c("ID_1 ID_2 missing sex\n0 0 0 D", paste0(data$S, " ", data$S, " ", 0, " ", 1)), fileConn)
      close(fileConn)

      # convert to BGEN file & plink
      system("qctool -g data/genotypes.gen -og data/genotypes.bgen")
      system("bgenix -g data/genotypes.bgen -clobber -index")
      system("qctool -g data/genotypes.gen -s data/samples.txt -og data/genotypes -ofiletype binary_ped")
      system("sed 's/^/S/g' -i data/genotypes.fam")
      
      # run vGWAS
      system("../build/bin/jlst_cpp -v data/phenotypes.csv -s , -o data/gwas_bgen.txt -g data/genotypes.bgen -p Y -i S -t 1")
      system("../build/bin/jlst_cpp -v data/phenotypes.csv -s , -o data/gwas_plink.txt -g data/genotypes -p Y -i S -t 1")
      system("osca --vqtl --bfile data/genotypes --pheno data/phenotypes.txt --out data/osca.txt --vqtl-mtd 1")

      # B-P using R implementation
      res_r <- bp(data$X, data$Y)

      # B-P using C++
      res_cpp_bgen <- fread("data/gwas_bgen.txt", select = c("beta", "se", "p", "phi_x", "se_x", "phi_xsq", "se_xsq", "phi_p"), col.names = c("BETA_mu.cpp_bgen", "SE_mu.cpp_bgen", "P_mu.cpp_bgen", "BETA_x.cpp_bgen", "SE_x.cpp_bgen", "BETA_xsq.cpp_bgen", "SE_xsq.cpp_bgen", "P.cpp_bgen"))
      res_cpp_plink <- fread("data/gwas_plink.txt", select = c("beta", "se", "p", "phi_x", "se_x", "phi_xsq", "se_xsq", "phi_p"), col.names = c("BETA_mu.cpp_plink", "SE_mu.cpp_plink", "P_mu.cpp_plink", "BETA_x.cpp_plink", "SE_x.cpp_plink", "BETA_xsq.cpp_plink", "SE_xsq.cpp_plink", "P.cpp_plink"))
      
      # Levene using OSCA
      res_osca <- fread("data/osca.txt.vqtl", select = c("beta", "se", "P"), col.names = c("BETA_x.osca", "SE_x.osca", "P.osca"))

      # combine results
      res <- cbind(res_r, res_cpp_bgen, res_cpp_plink, res_osca)

      # add LM
      fit <- tidy(lm(Y ~ X * U, data=data))
      res$BETA.x <- fit$estimate[2]
      res$BETA.u <- fit$estimate[3]
      res$BETA.xu <- fit$estimate[4]
      res$SE.x <- fit$std.error[2]
      res$SE.u <- fit$std.error[3]
      res$SE.xu <- fit$std.error[4]

      fit <- tidy(lm(Y ~ X, data=data))
      res$BETA_mu.r <- fit$estimate[2]
      res$SE_mu.r <- fit$std.error[2]
      res$P_mu.r <- fit$p.value[2]

      # add expected variance parameters
      res$EXP_x <- 2*delta*theta
      res$EXP_xsq <- theta*theta

      # add params
      res$phi <- phi
      res$af <- af
      res$lambda <- lambda
      res$theta <- theta
      res$delta <- delta

      # store result
      results <- rbind(results, res)
    }
  }

}

# save results for plotting
write.csv(results, file = paste0("data/results_", opt$dist, ".csv"))