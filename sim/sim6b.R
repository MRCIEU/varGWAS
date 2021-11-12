library("broom")
library("dplyr")
library("car")
source("funs.R")
set.seed(1243)

n_obs <- 10000
n_sim <- 300

# delta method does not give the correct SE for var(Y|G==2)
dummy_model_delta <- function(x, y, covar=NULL){
    if (!is.null(covar)){
        X <- as.matrix(cbind(x, covar))
    } else {
        X <- as.matrix(data.frame(x))
    }
    # first-stage fit
    fit <- qrfit(X=X, y=y, tau=.5, method="mm")
    b <- rbind(fit$b, fit$beta)
    # predicted
    X <- cbind(rep(1, nrow(X)), X)
    fitted <- X %*% b
    # residual
    d <- y - fitted
    # abs residual
    d <- abs(as.vector(d))
    # dummy SNP
    x <- as.factor(x)
    # second-stage model
    fit2 <- lm(d ~ x)
    # deltamethod
    v1 <- deltaMethod(fit2, "(2*b0*b1+b1^2)/(2/pi)", parameterNames=c("b0", "b1", "b2"))
    v2 <- deltaMethod(fit2, "(2*b0*b2+b2^2)/(2/pi)", parameterNames=c("b0", "b1", "b2"))
    return(
        c(v1$Estimate, v2$Estimate, v1$SE, v2$SE)
    )
}

df <- data.frame()
for (b in seq(0,6,2)){
    for (i in 1:n_sim){
        message(paste0("b:", b, " i:", i))
        x <- get_simulated_genotypes(0.4, n_obs)
        u <- rnorm(n_obs)
        y <- x*u*b + rnorm(n_obs)

        # estimate variance effect + SE
        res <- dummy_model_delta(x, y)

        df <- rbind(df, data.frame(
            b, t1=var(y[x==1]) - var(y[x==0]), t2=var(y[x==2]) - var(y[x==0]),
            b1=res[1], b2=res[2],
            s1=res[3], s2=res[4]
        ))
    }
}

df$lci1 <- df$b1 - (1.96 * df$s1)
df$lci2 <- df$b2 - (1.96 * df$s2)
df$uci1 <- df$b1 + (1.96 * df$s1)
df$uci2 <- df$b2 + (1.96 * df$s2)

df$t1_t <- 1^2*df$b^2
df$t2_t <- 2^2*df$b^2

r1 <- df %>% 
    dplyr::group_by(b) %>%
    dplyr::summarize(binom.test(sum(lci1 <= t1_t & uci1 >= t1_t), n()) %>% tidy)
r2 <- df %>% 
    dplyr::group_by(b) %>%
    dplyr::summarize(binom.test(sum(lci2 <= t2_t & uci2 >= t2_t), n()) %>% tidy)
r1 %>% dplyr::filter(conf.low > .95 | conf.high < .95)
r2 %>% dplyr::filter(conf.low > .95 | conf.high < .95)