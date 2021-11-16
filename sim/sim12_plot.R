library("ggplot2")
library("data.table")
set.seed(123)

# load in simulated results
results <- fread("results.csv")

# estimate mean and 95% CI
b1_dummy_ci <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarize(t.test(b1_dummy) %>% tidy) %>%
    dplyr::select(phi, estimate, conf.low, conf.high) %>%
    dplyr::mutate(genotype=1, method="LAD-BF")
b2_dummy_ci <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarize(t.test(b2_dummy) %>% tidy) %>%
    dplyr::select(phi, estimate, conf.low, conf.high) %>%
    dplyr::mutate(genotype=2, method="LAD-BF")
b1_osca_ci <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarize(t.test(b1_osca) %>% tidy) %>%
    dplyr::select(phi, estimate, conf.low, conf.high) %>%
    dplyr::mutate(genotype=1, method="OSCA")
b2_osca_ci <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarize(t.test(b2_osca) %>% tidy) %>%
    dplyr::select(phi, estimate, conf.low, conf.high) %>%
    dplyr::mutate(genotype=2, method="OSCA")
v1_mean <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarize(x=mean(v1-v0))
v2_mean <- results %>%
    dplyr::group_by(phi) %>%
    dplyr::summarize(x=mean(v2-v0))
v1 <- rbind(b1_dummy_ci, b1_osca_ci)
v1 <- merge(v1, v1_mean, "phi")
v2 <- rbind(b2_dummy_ci, b2_osca_ci)
v2 <- merge(v2, v2_mean, "phi")
ci <- rbind(v1, v2)
ci$phi <- as.factor(ci$phi)
ci$genotype <- as.factor(ci$genotype)

# plot
# TODO contrast with linear effect of X on var(Y)
pdf("sim12.pdf")
ggplot(data=ci, aes(x=x, y=estimate, ymin=conf.low, ymax=conf.high, group=genotype, shape=genotype)) +
    geom_point() + 
    geom_smooth(method="lm", formula=y~0+x, se=F, linetype = "dashed", color="grey") +
    geom_errorbar(width=.05) +
    theme_classic() +
    xlab("Difference in variance compared with G=0") +
    ylab("Estimated difference in variance compared with G=0 (95% CI)") +
    facet_grid(. ~ method) +
    labs(shape="Genotype") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    theme(
        strip.background = element_blank(),
        legend.title.align=0.5
    )
dev.off()