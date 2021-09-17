library("pwr")

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

get_osca <- function(x, y){
    s <- paste0("S", seq(1, length(y)))

    # write out GEN file
    write_gen("data/genotypes.gen", "01", "SNPID_1", "RSID_1", "1", "A", "G", x)

    # write phenotype & sample file
    write.table(file = "data/phenotypes.txt", sep = "\t", quote = F, row.names = F, col.names = F, data.frame(s, s, y))
    fileConn <- file("data/samples.txt")
    writeLines(c("ID_1 ID_2 missing sex\n0 0 0 D", paste0(s, " ", s, " ", 0, " ", 1)), fileConn)
    close(fileConn)

    # convert to BGEN file & plink
    system("qctool -g data/genotypes.gen -s data/samples.txt -og data/genotypes -ofiletype binary_ped", ignore.stdout = T,ignore.stderr = T)
    system("awk '{print \"S\"$0}' data/genotypes.fam > data/t; mv data/t data/genotypes.fam", ignore.stdout = T,ignore.stderr = T)

    # run OSCA
    system("osca --vqtl --bfile data/genotypes --pheno data/phenotypes.txt --out data/osca.txt --vqtl-mtd 1", ignore.stdout = T,ignore.stderr = T)
    osca <- fread("data/osca.txt.vqtl")

    return(osca)
}