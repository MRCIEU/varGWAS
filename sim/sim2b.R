library("dplyr")
library("broom")
library("tidyr")
library("ggpubr")
library("lmtest")
library("jlst")
library("data.table")
source("funs.R")
set.seed(123)

# Requires OSCA and QCTOOL on PATH

n_sim <- 1000
n_obs <- 100000
af <- 0.05
b <- 0

results <- data.frame()
for (dist in c("Normal", "T", "Lognormal", "Mixed Normal")){
    for (i in 1:n_sim){
        # simulate covariates
        data <- data.frame(
            S = paste0("S", seq(1, n_obs)),
            X = get_simulated_genotypes(af, n_obs),
            stringsAsFactors=F
        )

        if (dist == "Normal"){
            data$Y <- data$X * b + rnorm(n_obs)
        } else if (dist == "T"){
            data$Y <- data$X * b + rt(n_obs, 4)
        } else if (dist == "Lognormal"){
            data$Y <- data$X * b + rlnorm(n_obs)
        } else if (dist == "Mixed Normal"){
            data$Y <- data$X * b + c(rnorm(n_obs * .9), rnorm(n_obs * .1, mean=5))
        }

        # run models
        res <- run_models(data)
        bp_p <- vartest(data$Y, data$X, covar=NULL, covar.var=F, type=1, x.sq=T)$test$P
        res$bp_p <- bp_p
        res$dist <- dist
        res$af <- af
        res$b <- b

        # store result
        results <- rbind(results, res)
    }
}

qqgplot <- function(data, af, pcol, ci = 0.95) {
    temp <- data.frame()

    for (dist in c("Normal", "T", "Lognormal", "Mixed Normal")){
        p <- data %>%
            dplyr::filter(dist == !!dist & af == !!af) %>%
            tidyr::drop_na(!!pcol) %>%
            dplyr::pull(!!pcol)
        n  <- length(p)
        temp <- rbind(temp, data.frame(
            dist,
            observed = -log10(sort(p)),
            expected = -log10(ppoints(n)),
            clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
            cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
        ))
    }

    log10Pe <- expression(paste("Expected -log"[10], plain(P)))
    log10Po <- expression(paste("Observed -log"[10], plain(P)))

    pl <- ggplot(temp) +
        geom_point(aes(expected, observed), shape = 1, size = 1) +
        geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
        geom_line(aes(expected, cupper), linetype = 2) +
        geom_line(aes(expected, clower), linetype = 2) +
        xlab(log10Pe) +
        ylab(log10Po) +
        facet_wrap(vars(dist), scale="free_y") +
        theme_minimal()

    return(pl)
}

# print warnings
warnings()

p1 <- qqgplot(results, 0.05, "P.osca_median")
p2 <- qqgplot(results, 0.05, "bp_p")
p3 <- qqgplot(results, 0.05, "P.cpp_bf")

p <- ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol = 3, nrow = 1)
pdf("data/t1e_10k.pdf", height=7, width=21)
print(p)
dev.off()