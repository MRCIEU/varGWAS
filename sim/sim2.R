library("dplyr")
library("varGWASR")
library("data.table")
library("jlst")
library("ggplot2")
library("ggpubr")
source("funs.R")
set.seed(123)

n_sim <- 1000
n_obs <- 10000
b <- 0

results <- data.frame()
for (af in c(0.05, 0.1, 0.15, 0.2)){
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

            bp_p <- vartest(data$Y, data$X, covar=NULL, covar.var=F, type=1, x.sq=T)$test$P
            bf_p <- varGWASR::model(data, "X", "Y", covar1 = NULL, covar2 = NULL)$phi_p
            res <- data.frame(
                bp_p, bf_p
            )
            res$dist <- dist
            res$af <- af
            res$b <- b

            # store result
            results <- rbind(results, res)
        }
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

# save results for plotting
write.csv(results, file = paste0("data/sim2.csv"))

p1 <- qqgplot(results, 0.05, "bp_p")
p2 <- qqgplot(results, 0.1, "bp_p")
p3 <- qqgplot(results, 0.15, "bp_p")
p4 <- qqgplot(results, 0.2, "bp_p")

p <- ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
pdf("data/maf_t1e_bp.pdf")
print(p)
dev.off()

p1 <- qqgplot(results, 0.05, "bf_p")
p2 <- qqgplot(results, 0.1, "bf_p")
p3 <- qqgplot(results, 0.15, "bf_p")
p4 <- qqgplot(results, 0.2, "bf_p")

p <- ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
pdf("data/maf_t1e_bf.pdf")
print(p)
dev.off()