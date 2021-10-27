library("car")
library("quantreg")
set.seed(124)

n_sim <- 1000
n_obs <- 1000

#' vartest
#'
#' vartest performs variability tests by either the Breusch-Pagan or Brown-Forsythe methods.
#' @param y vector of outcome values.
#' @param x vector of exposure values.
#' @param covar a data.frame of covariates.
#' @param covar.var adjust the second stage (variance component) of the approach by the covariates.
#' @param type type of test (default: 1 [Breusch-Pagan variance test]; options: 1 [Breusch-Pagan variance test], 2 [Brown-Forsythe variance test]).
#' @param x.sq include x-squared in the variance part of the model.
#' @return a list of results. F is the test statistic, DF is the degrees of freedom and P is the p-value. The model coefficients from variance part of the model are given in the coef object.
#' @examples
#' x <- rbinom(1000, 1, 0.5)
#' y <- 0.5 + 0.025*x + rnorm(1000, 0, sqrt(0.005*x)) + rnorm(1000, 0, 0.1)
#' vartest(y, x, type=2)
#' @author James R Staley <jrstaley95@gmail.com>
#' @export
vartest <- function(y, x, covar=NULL, covar.var=F, type=1, x.sq=F){
  
  # Errors
  if(!(is.numeric(y) | is.integer(y))) stop("y has to be a numeric variable")
  if(!(is.numeric(x) | is.integer(x) | is.factor(x))) stop("x has to be either a numeric variable or a factor variable")
  if(!is.null(covar)){if(!is.data.frame(covar)) stop("covar has to be a data.frame")}
  if(is.null(covar) & covar.var) stop("covar.var cannot be TRUE if there are no covariates")
  if(length(y)!=length(x)) stop("y is not the same size as x")
  if(!is.null(covar)){if(length(y)!=nrow(covar)) stop("y is not the same size as covar")}
  if(!(type %in% 1:2)) stop("type has to be set to either 1 or 2")
  if(is.logical(x.sq)==F) stop("x.sq has to logical")
  if(x.sq==T & is.factor(x)==T) stop("x.sq cannot be set to true if x is a factor variable")
  
  # Missing values
  data <- cbind(y, x); if(!is.null(covar)){data <- cbind(data, covar)}
  keep <- complete.cases(data)
  y <- y[keep]; x <- x[keep]; if(!is.null(covar)){covar <- covar[keep,,drop=F]}
  
  # Covariates
  if(!is.null(covar)){
    covar <- model.matrix(as.formula(~ .), data=covar)[,-1,drop=F]
    if(any(is.na(covar))) stop("there are missing values in the covariates")
  }
  
  # Exposure squared
  if(x.sq==T){
    x2 <- x^2
  }
  
  # Variance test
  if(!is.null(covar)){
    if(type==1){d <- (abs(resid(lm(y~x+covar))))}else{d <- suppressWarnings(abs(resid(rq(y~x+covar, tau=0.5))))}
  }else{
    if(type==1){d <- (abs(resid(lm(y~x))))}else{d <- suppressWarnings(abs(resid(rq(y~x, tau=0.5))))}
  }
  if(!is.null(covar) & covar.var){
    if(x.sq==T){mod <- lm(d~x+x2+covar); mod0 <- lm(d~covar)}else{mod <- lm(d~x+covar); mod0 <- lm(d~covar)}
  }else{
    if(x.sq==T){mod <- lm(d~x+x2); mod0 <- lm(d~1)}else{mod <- lm(d~x); mod0 <- lm(d~1)}
  }
  coef <- summary(mod)$coefficients; rownames(coef) <- sub("covar", "", rownames(coef))
  test <- anova(mod0, mod)[2,c(3,5,6)]; test <- as.data.frame(test); names(test) <- c("DF", "F", "P"); rownames(test) <- 1:nrow(test)
  results <- list(coef=coef, test=test)
  
  # Results  
  return(results)
  
}

p_levene <- rep(NA, n_sim)
p_bp <- rep(NA, n_sim)
p_glejser <- rep(NA, n_sim)
for (i in 1:n_sim){
    x <- rbinom(n_obs, 2, .5)
    u <- rnorm(n_obs)
    y <- x*u*.5 + rnorm(n_obs)
    p_bp[i] <- vartest(y, x, type = 1, x.sq=T)$test$P
    p_levene[i] <- leveneTest(y, x, center = mean)$P
    p_glejser[i] <- glejser(lm(y ~ x))$p.value
}