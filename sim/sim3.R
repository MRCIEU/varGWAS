library("dplyr")
library("broom")
library("tidyr")
library("ggpubr")
source("funs.R")
set.seed(123)

n_sim <- 1000
n_obs <- 1000

results <- data.frame()
for (trans in c("log", "sqrt", "irnt", "cube_root")){
    for (dist in c("Normal", "T", "Lognormal", "Mixed Normal")){
        for (i in 1:n_sim){
            # effect size for SNP to have 80% power
            b <- 0.2125
            x <- get_simulated_genotypes(0.1, n_obs)
            if (dist == "Normal"){
                y <- 100 + x * b + rnorm(n_obs)
            } else if (dist == "T"){
                y <- 100 + x * b + rt(n_obs, 4)
            } else if (dist == "Lognormal"){
                y <- 100 + x * b + rlnorm(n_obs)
            } else if (dist == "Mixed Normal"){
                y <- 100 + x * b + c(rnorm(n_obs * .9), rnorm(n_obs * .1, mean=5))
            }

            # test for main effect
            lm_p <- tidy(lm(y ~ x))$p.value[2]

            if (trans == "log"){
                y <- log(y)
            } else if (trans == "sqrt"){
                y <- sqrt(y)
            } else if (trans == "irnt"){
                y <- irnt(y)
            } else if (trans == "cube_root"){
                y <- y^(1/3)
            }

            # test for effect using B-P
            res <- data.frame(
                dist,
                trans,
                lm_p,
                b,
                bp_p=vartest(y, x, type=1, x.sq=T)$test$P,
                osca_p=get_osca(x, y)$P
            )
            results <- rbind(results, res)
        }
    }
}

# check main effect has 80% power
results %>% group_by(dist) %>% summarize(tidy(binom.test(sum(lm_p < 0.05), n())))

qqgplot <- function(data, trans, pcol, ci = 0.95) {
    temp <- data.frame()

    for (dist in c("Normal", "T", "Lognormal", "Mixed Normal")){
        p <- data %>%
            filter(dist == !!dist & trans == !!trans) %>%
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

p1 <- qqgplot(results, "log", "bp_p")
p2 <- qqgplot(results, "sqrt", "bp_p")
p3 <- qqgplot(results, "irnt", "bp_p")
p4 <- qqgplot(results, "cube_root", "bp_p")

p <- ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
pdf("data/trans_t1e_bp.pdf")
print(p)
dev.off()

p1 <- qqgplot(results, "log", "osca_p")
p2 <- qqgplot(results, "sqrt", "osca_p")
p3 <- qqgplot(results, "irnt", "osca_p")
p4 <- qqgplot(results, "cube_root", "osca_p")

p <- ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
pdf("data/trans_t1e_osca.pdf")
print(p)
dev.off()
