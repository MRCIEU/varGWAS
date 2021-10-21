library("pwr")
library("jlst")

#' Function to simulate genotypes in HWE
#' @param q Recessive/alternative allele frequency
#' @param n_obs Number of observations to return
get_simulated_genotypes <- function(q, n_obs){
  p <- 1 - q
  x <- sample(c(0, 1, 2), n_obs, prob=c(p^2, 2 * p * q, q^2), replace=T)
  return(x)
}

#' Estimate treatment effect size a binary exposure
#' @param n_obs Number of samples in analysis
#' @param treatment1_p The probability of recieving treatment
#' @param sd The SD out of the normal outcome
#' @param alpha P value threshold
#' @param power Power to detect effect
#' @return Effect size
get_binary_delta <- function(n_obs, treatment1_p, sd, alpha, power){
  p <- power.t.test(n=n_obs * treatment1_p, delta=NULL, sd=sd, sig.level=alpha, power=power, type = c("two.sample"), alternative = c("two.sided"))
  return(p$delta)
}

#' Estimate treatment effect size given a continuous exposure
#' @param n_obs Number of samples in analysis
#' @param alpha P value threshold
#' @param power Power to detect effect
#' @return Effect size
get_cont_delta <- function(n_obs, alpha, power){
  # u = number of terms in the model (excluding intercept)
  u <- 1

  # v = error degrees of freedom
  v <- n_obs - u - 1

  # f2 = Cohen f2 (variance explained by the model)
  p <- pwr.f2.test(u = u, v = v, sig.level = alpha, power = power)

  return(sqrt(p$f2))
}

write_gen <- function(f, chr, snpid, rsid, pos, A1, A2, x){
  fileConn <- file(f)
  writeLines(c(paste(chr, snpid, rsid, pos, A1, A2, paste(sapply(x, function(g) if (g == 0) { "0 0 1" } else if (g == 1) { "0 1 0" } else if (g == 2) { "1 0 0" }), collapse = " "), collapse = " ")), fileConn)
  close(fileConn)
}

run_models <- function(data){
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
  bp.time <- system.time(system("varGWAS -v data/phenotypes.csv -s , -o data/gwas-bp.txt -b data/genotypes.bgen -p Y -i S -t 1"))
  bf.time <- system.time(system("varGWAS -v data/phenotypes.csv -s , -o data/gwas-bf.txt -b data/genotypes.bgen -p Y -i S -t 1 -r"))
  osca_mean.time <- system.time(system("osca --vqtl --bfile data/genotypes --pheno data/phenotypes.txt --out data/osca-mean.txt --vqtl-mtd 1"))
  osca_median.time <- system.time(system("osca --vqtl --bfile data/genotypes --pheno data/phenotypes.txt --out data/osca-median.txt --vqtl-mtd 2"))
  bp.time <- t(data.matrix(bp.time)) %>% as.data.frame
  bf.time <- t(data.matrix(bf.time)) %>% as.data.frame
  osca_mean.time <- t(data.matrix(osca_mean.time)) %>% as.data.frame
  osca_median.time <- t(data.matrix(osca_median.time)) %>% as.data.frame
  names(bp.time) <- paste0(names(bp.time), ".cpp_bp")
  names(bf.time) <- paste0(names(bf.time), ".cpp_bf")
  names(osca_mean.time) <- paste0(names(osca_mean.time), ".osca_mean")
  names(osca_median.time) <- paste0(names(osca_median.time), ".osca_median")

  # R
  bp <- vartest(data$Y, x = data$X, type = 1, x.sq = T)
  bf <- vartest(data$Y, x = data$X, type = 2, x.sq = T)
  bp$coef <- rbind(bp$coef, c(NA, NA, NA, NA))
  bf$coef <- rbind(bf$coef, c(NA, NA, NA, NA))
  res_r_bp <- data.frame(BETA_x.r_bp = bp$coef[2, 1], SE_x.r_bp = bp$coef[2, 2], BETA_xsq.r_bp = bp$coef[3, 1], SE_xsq.r_bp = bp$coef[3, 2], P.r_bp = as.numeric(bp$test[3]))
  res_r_bf <- data.frame(BETA_x.r_bf = bf$coef[2, 1], SE_x.r_bf = bf$coef[2, 2], BETA_xsq.r_bf = bf$coef[3, 1], SE_xsq.r_bf = bf$coef[3, 2], P.r_bf = as.numeric(bf$test[3]))

  # C++
  res_cpp_bp <- fread("data/gwas-bp.txt", select = c("beta", "beta_lad", "se", "p", "phi_x", "se_x", "phi_xsq", "se_xsq", "phi_p"), col.names = c("BETA_mu.cpp_bp", "BETA_lad.cpp_bp", "SE_mu.cpp_bp", "P_mu.cpp_bp", "BETA_x.cpp_bp", "SE_x.cpp_bp", "BETA_xsq.cpp_bp", "SE_xsq.cpp_bp", "P.cpp_bp"))
  res_cpp_bf <- fread("data/gwas-bf.txt", select = c("beta", "beta_lad", "se", "p", "phi_x", "se_x", "phi_xsq", "se_xsq", "phi_p"), col.names = c("BETA_mu.cpp_bf", "BETA_lad.cpp_bf", "SE_mu.cpp_bf", "P_mu.cpp_bf", "BETA_x.cpp_bf", "SE_x.cpp_bf", "BETA_xsq.cpp_bf", "SE_xsq.cpp_bf", "P.cpp_bf"))

  if (nrow(res_cpp_bp) == 0){
    res_cpp_bp <- data.frame(
      BETA_mu.cpp_bp=NA, BETA_lad.cpp_bp=NA, SE_mu.cpp_bp=NA, P_mu.cpp_bp=NA, BETA_x.cpp_bp=NA, SE_x.cpp_bp=NA, BETA_xsq.cpp_bp=NA, SE_xsq.cpp_bp=NA, P.cpp_bp=NA
    )
  }

  if (nrow(res_cpp_bf) == 0){
    res_cpp_bf <- data.frame(
      BETA_mu.cpp_bf=NA, BETA_lad.cpp_bf=NA, SE_mu.cpp_bf=NA, P_mu.cpp_bf=NA, BETA_x.cpp_bf=NA, SE_x.cpp_bf=NA, BETA_xsq.cpp_bf=NA, SE_xsq.cpp_bf=NA, P.cpp_bf=NA
    )
  }

  # Levene using OSCA
  # Note this method will not produce effect or se if P==0
  res_osca_mean <- fread("data/osca-mean.txt.vqtl", select = c("beta", "se", "P"), col.names = c("BETA_x.osca_mean", "SE_x.osca_mean", "P.osca_mean"))
  res_osca_median <- fread("data/osca-median.txt.vqtl", select = c("beta", "se", "P"), col.names = c("BETA_x.osca_median", "SE_x.osca_median", "P.osca_median"))

  # combine results
  res <- cbind(res_r_bp, res_r_bf, res_cpp_bp, res_cpp_bf, res_osca_mean, res_osca_median, bp.time, bf.time, osca_mean.time, osca_median.time)

  # add LM
  if ("U" %in% names(data)){
    fit <- tidy(lm(Y ~ X * U, data = data))
    res$BETA.x <- fit$estimate[2]
    res$BETA.u <- fit$estimate[3]
    res$BETA.xu <- fit$estimate[4]
    res$SE.x <- fit$std.error[2]
    res$SE.u <- fit$std.error[3]
    res$SE.xu <- fit$std.error[4]
  }

  fit <- tidy(lm(Y ~ X, data = data))
  res$BETA_mu.r <- fit$estimate[2]
  res$SE_mu.r <- fit$std.error[2]
  res$P_mu.r <- fit$p.value[2]
  res$n <- nrow(data)
  
  return(res)
}

# taken from https://github.com/MRCIEU/PHESANT/blob/3f4a65d7fe93aaf01f3a4a3f39843562612a8d65/WAS/testContinuous.r#L243
irnt <- function(pheno) {
	numPhenos = length(which(!is.na(pheno)))
	quantilePheno = (rank(pheno, na.last="keep", ties.method="random")-0.5)/numPhenos
	phenoIRNT = qnorm(quantilePheno)	
	return(phenoIRNT);
}