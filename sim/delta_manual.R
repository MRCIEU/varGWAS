library("car")
set.seed(123)

n_obs <- 1000000
n_sim <- 200

get_se <- function(grad, fit){
    vb <- vcov(fit)
    vG <- t(grad) %*% vb %*% grad
    sqrt(vG)
}

results <- data.frame()
for (i in 1:n_sim){
    b0 <- rnorm(1)
    b1 <- rnorm(1)
    b2 <- rnorm(1)
    x <- rbinom(n_obs, 2, runif(1))
    x1 <- as.numeric(x==1)
    x2 <- as.numeric(x==2)
    d <- b0 + x1*b1 + x2*b2 + rnorm(n_obs)
    fit2 <- lm(d ~ x1 + x2)

    # delta
    tryCatch({
        se_b0 <- get_se(c(
            2 * b0/(2/pi),
            0,
            0
        ), fit2)

        se_b1 <- get_se(c(
            2 * b0/(2/pi),
            (2 * b0 + 2 * b1)/(2/pi),
            0
        ), fit2)

        se_b2 <- get_se(c(
            2 * b0/(2/pi),
            0,
            (2 * b0 + 2 * b2)/(2/pi)
        ), fit2)
        
        delta_b0 <- deltaMethod(fit2, "b0^2/(2/pi)", parameterNames=c("b0", "b1", "b2"))
        delta_b1 <- deltaMethod(fit2, "b0^2/(2/pi) + (2*b0*b1+b1^2)/(2/pi)", parameterNames=c("b0", "b1", "b2"))
        delta_b2 <- deltaMethod(fit2, "b0^2/(2/pi) + (2*b0*b2+b2^2)/(2/pi)", parameterNames=c("b0", "b1", "b2"))

        results <- rbind(results, data.frame(
            est_se_b0=se_b0,
            est_se_b1=se_b1,
            est_se_b2=se_b2,
            est_beta_b0=b0^2/(2/pi),
            est_beta_b1=b0^2/(2/pi) + (2*b0*b1+b1^2)/(2/pi),
            est_beta_b2=b0^2/(2/pi) + (2*b0*b2+b2^2)/(2/pi),
            delta_se_b0=delta_b0$SE,
            delta_se_b1=delta_b1$SE,
            delta_se_b2=delta_b2$SE,
            delta_beta_b0=delta_b0$Estimate,
            delta_beta_b1=delta_b1$Estimate,
            delta_beta_b2=delta_b2$Estimate
        ))

    }, error=function(e){

    })

}