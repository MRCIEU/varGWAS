library("ggplot2")
library("data.table")
library("dplyr")
library("broom")
set.seed(123)

# load in simulated results
results <- fread("results.csv")

# estimate mean and 95% CI
b1_dummy_ci <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarize(t.test(b1_dummy) %>% tidy) %>%
    dplyr::select(phi, estimate, conf.low, conf.high) %>%
    dplyr::mutate(genotype="SNP=1", method="LAD-BF")
b2_dummy_ci <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarize(t.test(b2_dummy) %>% tidy) %>%
    dplyr::select(phi, estimate, conf.low, conf.high) %>%
    dplyr::mutate(genotype="SNP=2", method="LAD-BF")
v1_mean <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarize(x=mean(v1))
v2_mean <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarize(x=mean(v2))
v1 <- b1_dummy_ci
v1 <- merge(v1, v1_mean, "phi")
v2 <- b2_dummy_ci
v2 <- merge(v2, v2_mean, "phi")
ci <- rbind(v1, v2)
ci$phi <- as.factor(ci$phi)
ci$genotype <- as.factor(ci$genotype)

# plot
pdf("sim12_estimate.pdf")
ggplot(data=ci, aes(x=x, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() + 
    geom_abline(intercept = 0, slope = 1, linetype="dashed", color="grey") +
    geom_errorbar(width=.05) +
    theme_classic() +
    xlab("Difference in variance compared with G=0") +
    ylab("Estimated difference in variance compared with G=0 (95% CI)") +
    facet_wrap(method ~ genotype, scales="free") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    theme(
        strip.background = element_blank(),
        legend.title.align=0.5
    )
dev.off()

# coverage
results <- merge(results, v1_mean, "phi")
results <- merge(results, v2_mean, "phi")
names(results)[10] <- "v1_mean"
names(results)[11] <- "v2_mean"
results$b1_dummy_lci <- results$b1_dummy - (1.96 * results$s1_dummy)
results$b1_dummy_uci <- results$b1_dummy + (1.96 * results$s1_dummy)
results$b2_dummy_lci <- results$b2_dummy - (1.96 * results$s2_dummy)
results$b2_dummy_uci <- results$b2_dummy + (1.96 * results$s2_dummy)

b1_dummy <- results %>% 
    dplyr::group_by(phi) %>%
    dplyr::summarize(binom.test(sum(b1_dummy_lci <= v1_mean & b1_dummy_uci >= v1_mean), n()) %>% tidy) %>%
    dplyr::select(phi, estimate, conf.low, conf.high) %>%
    dplyr::mutate(genotype="SNP=1", method="LAD-BF")
b2_dummy <- results %>% 
    dplyr::group_by(phi) %>%
    dplyr::summarize(binom.test(sum(b2_dummy_lci <= v2_mean & b2_dummy_uci >= v2_mean), n()) %>% tidy) %>%
    dplyr::select(phi, estimate, conf.low, conf.high) %>%
    dplyr::mutate(genotype="SNP=2", method="LAD-BF")

coverage1 <- b1_dummy
coverage1 <- merge(coverage1, v1_mean, "phi")
coverage2 <- b2_dummy
coverage2 <- merge(coverage2, v2_mean, "phi")
coverage <- rbind(coverage1, coverage2)
coverage$phi <- as.factor(coverage$phi)
coverage$genotype <- as.factor(coverage$genotype)

pdf("sim12_coverage.pdf")
ggplot(data=coverage, aes(x=x, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() + 
    geom_errorbar(width=.05) +
    theme_classic() +
    xlab("Difference in variance compared with G=0") +
    ylab("Coverage of 95% CI (95% CI)") +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey") +
    facet_grid(method ~ genotype, scales="free_x") +
    labs(shape="Genotype") +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks(n = 5)) +
    theme(
        strip.background = element_blank(),
        legend.title.align=0.5
    )
dev.off()