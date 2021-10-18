library("dplyr")
library("broom")
library("tidyr")
library("data.table")
library("ggpubr")
source("funs.R")
set.seed(123)

# Requires OSCA and QCTOOL on PATH

n_sim <- 1000
n_obs <- 1000

# effect size for SNP to have 80% power
b <- 0.2125

results <- data.frame()
for (trans in c("log", "sqrt", "irnt", "cube_root")){
    for (dist in c("Normal", "T", "Lognormal", "Mixed Normal")){
        for (i in 1:n_sim){
            # simulate covariates
            data <- data.frame(
                S = paste0("S", seq(1, n_obs)),
                X = get_simulated_genotypes(0.1, n_obs),
                stringsAsFactors=F
            )

            # simulate outcome
            data$Y <- 100 + data$X * b
            if (dist == "Normal"){
                data$Y <- data$Y + rnorm(n_obs)
            } else if (dist == "T"){
                data$Y <- data$Y + rt(n_obs, 4)
            } else if (dist == "Lognormal"){
                data$Y <- data$Y + rlnorm(n_obs)
            } else if (dist == "Mixed Normal"){
                data$Y <- data$Y + c(rnorm(n_obs * .9), rnorm(n_obs * .1, mean=5))
            }

            # test for main effect
            lm_p <- tidy(lm(Y ~ X, data=data))$p.value[2]

            if (trans == "log"){
                data$Y <- log(data$Y)
            } else if (trans == "sqrt"){
                data$Y <- sqrt(data$Y)
            } else if (trans == "irnt"){
                data$Y <- irnt(data$Y)
            } else if (trans == "cube_root"){
                data$Y <- data$Y^(1/3)
            }

            # run models
            res <- run_models(data)
            res$dist <- dist
            res$trans <- trans
            res$lm_p <- lm_p
            res$b <- b

            # store result
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
            drop_na(!!pcol) %>%
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

p1 <- qqgplot(results, "log", "P.cpp_bp")
p2 <- qqgplot(results, "sqrt", "P.cpp_bp")
p3 <- qqgplot(results, "irnt", "P.cpp_bp")
p4 <- qqgplot(results, "cube_root", "P.cpp_bp")

p <- ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
pdf("data/trans_t1e_bp.pdf")
print(p)
dev.off()

p1 <- qqgplot(results, "log", "P.cpp_bf")
p2 <- qqgplot(results, "sqrt", "P.cpp_bf")
p3 <- qqgplot(results, "irnt", "P.cpp_bf")
p4 <- qqgplot(results, "cube_root", "P.cpp_bf")

p <- ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
pdf("data/trans_t1e_bf.pdf")
print(p)
dev.off()

p1 <- qqgplot(results, "log", "P.osca")
p2 <- qqgplot(results, "sqrt", "P.osca")
p3 <- qqgplot(results, "irnt", "P.osca")
p4 <- qqgplot(results, "cube_root", "P.osca")

p <- ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
pdf("data/trans_t1e_osca.pdf")
print(p)
dev.off()