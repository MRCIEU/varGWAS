library("dplyr")
library("broom")
library("tidyr")
library("ggpubr")
library("lmtest")
library("data.table")
source("funs.R")
set.seed(123)

# Requires OSCA and QCTOOL on PATH

n_sim <- 1000
n_obs <- 1000
b <- 0

results <- data.frame()
for (af in c(0.01, 0.05, 0.1, 0.2)){
    for (dist in c("Normal", "T", "Lognormal", "Mixed Normal")){
        for (i in 1:n_sim){
            x <- get_simulated_genotypes(af, n_obs)
            x2 <- x^2

            if (dist == "Normal"){
                y <- x * b + rnorm(n_obs)
            } else if (dist == "T"){
                y <- x * b + rt(n_obs, 4)
            } else if (dist == "Lognormal"){
                y <- x * b + rlnorm(n_obs)
            } else if (dist == "Mixed Normal"){
                y <- x * b + c(rnorm(n_obs * .9), rnorm(n_obs * .1, mean=5))
            }

            # test for variance effect
            res <- data.frame(
                dist,
                af,
                b,
                bp_p=vartest(y, x, type=1, x.sq=T)$test$P,
                osca_p=get_osca(x, y)$P
            )
            results <- rbind(results, res)
        }
    }
}

qqgplot <- function(data, af, pcol, ci = 0.95) {
    temp <- data.frame()

    for (dist in c("Normal", "T", "Lognormal", "Mixed Normal")){
        p <- data %>%
            filter(dist == !!dist & af == !!af) %>%
            pull(!!pcol)
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

p1 <- qqgplot(results, 0.01, "bp_p")
p2 <- qqgplot(results, 0.05, "bp_p")
p3 <- qqgplot(results, 0.1, "bp_p")
p4 <- qqgplot(results, 0.2, "bp_p")

p <- ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
pdf("data/maf_t1e_bp.pdf")
print(p)
dev.off()

p1 <- qqgplot(results, 0.01, "osca_p")
p2 <- qqgplot(results, 0.05, "osca_p")
p3 <- qqgplot(results, 0.1, "osca_p")
p4 <- qqgplot(results, 0.2, "osca_p")

p <- ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
pdf("data/maf_t1e_osca.pdf")
print(p)
dev.off()