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
get_simulated_genotypes <- function(q, n_obs){
  p <- 1 - q
  x <- sample(c(0, 1, 2), n_obs, prob=c(p^2, 2 * p * q, q^2), replace=T)
  return(x)
}