library("data.table")
library("broom")
library("genpwr")
library('optparse')
library("jlst")
source("funs.R")
set.seed(12345)

option_list <- list(
  make_option(c("-d", "--dist"), type = "character", default = "n", help = "Outcome distribution", metavar = "character")
);
opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

n_obs <- 200
n_sim <- 200
af <- 0.4

# main effect size of X on Y detectable with 80% power
delta <- as.numeric(genpwr.calc(calc = "es", model = "linear", ge.interaction = NULL, N = n_obs, k = NULL, MAF = af, Power = 0.8, Alpha = 0.05, sd_y = 1, True.Model = "Additive", Test.Model = "Additive")$ES_at_Alpha_0.05)

# simulate GxE interaction effects and estimate power
results <- data.frame()
for (phi in seq(0, 6, 0.5)) {
  theta <- delta * phi

  for (lambda in c(1, 10, 100, 1000)) {
    for (i in 1:n_sim) {

      # simulate covariates
      data <- data.frame(
        S = paste0("S", seq(1, n_obs * lambda)),
        X = get_simulated_genotypes(af, n_obs * lambda),
        U = rnorm(n_obs * lambda)
      )

      # simulate outcome
      data$Y <- data$X * delta +
        data$U * delta +
        data$X * data$U * theta

      # add error term
      if (opt$dist == "n") {
        data$Y <- data$Y + rnorm(n_obs * lambda)
      } else if (opt$dist == "t") {
        data$Y <- data$Y + rt(n_obs * lambda, 4)
      } else if (opt$dist == "l") {
        data$Y <- data$Y + rlnorm(n_obs * lambda)
      } else {
        stop(paste0("Distribution not implemented: ", opt$dist))
      }

      # write out GEN file
      write_gen("data/genotypes.gen", "01", "SNPID_1", "RSID_1", "1", "A", "G", data$X)

      # write phenotype & sample file
      write.table(file = "data/phenotypes.csv", sep = ",", quote = F, row.names = F, data)
      write.table(file = "data/phenotypes.txt", sep = "\t", quote = F, row.names = F, col.names = F, data[, c("S", "S", "Y")])
      fileConn <- file("data/samples.txt")
      writeLines(c("ID_1 ID_2 missing sex\n0 0 0 D", paste0(data$S, " ", data$S, " ", 0, " ", 1)), fileConn)
      close(fileConn)

      # convert to BGEN file & plink
      system("qctool -g data/genotypes.gen -og data/genotypes.bgen")
      system("bgenix -g data/genotypes.bgen -clobber -index")
      system("qctool -g data/genotypes.gen -s data/samples.txt -og data/genotypes -ofiletype binary_ped")
      system("sed 's/^/S/g' -i data/genotypes.fam")

      # run vGWAS
      system("varGWAS -v data/phenotypes.csv -s , -o data/gwas-bp.txt -b data/genotypes.bgen -p Y -i S -t 1")
      system("varGWAS -v data/phenotypes.csv -s , -o data/gwas-bf.txt -b data/genotypes.bgen -p Y -i S -t 1 -r")
      system("osca --vqtl --bfile data/genotypes --pheno data/phenotypes.txt --out data/osca.txt --vqtl-mtd 1")

      # R
      bp <- vartest(data$Y, x = data$X, type = 1, x.sq = T)
      bf <- vartest(data$Y, x = data$X, type = 2, x.sq = T)
      res_r_bp <- data.frame(BETA_x.r_bp = bp$coef[2, 1], SE_x.r_bp = bp$coef[2, 2], BETA_xsq.r_bp = bp$coef[3, 1], SE_xsq.r_bp = bp$coef[3, 2], P.r_bp = as.numeric(bp$test[3]))
      res_r_bf <- data.frame(BETA_x.r_bf = bf$coef[2, 1], SE_x.r_bf = bf$coef[2, 2], BETA_xsq.r_bf = bf$coef[3, 1], SE_xsq.r_bf = bf$coef[3, 2], P.r_bf = as.numeric(bf$test[3]))

      # C++
      res_cpp_bp <- fread("data/gwas-bp.txt", select = c("beta", "se", "p", "phi_x", "se_x", "phi_xsq", "se_xsq", "phi_p"), col.names = c("BETA_mu.cpp_bp", "SE_mu.cpp_bp", "P_mu.cpp_bp", "BETA_x.cpp_bp", "SE_x.cpp_bp", "BETA_xsq.cpp_bp", "SE_xsq.cpp_bp", "P.cpp_bp"))
      res_cpp_bf <- fread("data/gwas-bf.txt", select = c("beta", "se", "p", "phi_x", "se_x", "phi_xsq", "se_xsq", "phi_p"), col.names = c("BETA_mu.cpp_bf", "SE_mu.cpp_bf", "P_mu.cpp_bf", "BETA_x.cpp_bf", "SE_x.cpp_bf", "BETA_xsq.cpp_bf", "SE_xsq.cpp_bf", "P.cpp_bf"))

      # Levene using OSCA
      # Note this method will not produce effect or se if P==0
      res_osca <- fread("data/osca.txt.vqtl", select = c("beta", "se", "P"), col.names = c("BETA_x.osca", "SE_x.osca", "P.osca"))

      # combine results
      res <- cbind(res_r_bp, res_r_bf, res_cpp_bp, res_cpp_bf, res_osca)

      # add LM
      fit <- tidy(lm(Y ~ X * U, data = data))
      res$BETA.x <- fit$estimate[2]
      res$BETA.u <- fit$estimate[3]
      res$BETA.xu <- fit$estimate[4]
      res$SE.x <- fit$std.error[2]
      res$SE.u <- fit$std.error[3]
      res$SE.xu <- fit$std.error[4]
      fit <- tidy(lm(Y ~ X, data = data))
      res$BETA_mu.r <- fit$estimate[2]
      res$SE_mu.r <- fit$std.error[2]
      res$P_mu.r <- fit$p.value[2]

      # add expected variance parameters
      res$EXP_x <- 2 * delta * theta
      res$EXP_xsq <- theta * theta

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

# print warnings
warnings()

# save results for plotting
write.csv(results, file = paste0("data/power_", opt$dist, ".csv"))
