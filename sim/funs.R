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