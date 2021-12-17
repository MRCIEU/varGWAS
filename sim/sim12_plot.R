library("ggplot2")
library("data.table")
library("dplyr")
library("broom")
library("ggpubr")
set.seed(123)

get_est <- function(file){
    # load in simulated results
    results <- fread(file)

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
    b1_osca_ci <- results %>%
        dplyr::group_by(phi) %>%
        dplyr::summarize(t.test(b1_osca) %>% tidy) %>%
        dplyr::select(phi, estimate, conf.low, conf.high) %>%
        dplyr::mutate(genotype="SNP=1", method="OSCA-BF")
    b2_osca_ci <- results %>%
        dplyr::group_by(phi) %>%
        dplyr::summarize(t.test(b2_osca) %>% tidy) %>%
        dplyr::select(phi, estimate, conf.low, conf.high) %>%
        dplyr::mutate(genotype="SNP=2", method="OSCA-BF")
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

    return(ci)
}

get_coverage <- function(file){
    # coverage
    results <- fread(file)

    # estimate differences in variance
    v1_mean <- results %>%
        dplyr::group_by(phi) %>%
        dplyr::summarize(x=mean(v1-v0))
    v2_mean <- results %>%
        dplyr::group_by(phi) %>%
        dplyr::summarize(x=mean(v2-v0))
    
    # coverage
    results <- merge(results, v1_mean, "phi")
    results <- merge(results, v2_mean, "phi")
    names(results)[15] <- "v1_mean"
    names(results)[16] <- "v2_mean"
    results$b1_dummy_lci <- results$b1_dummy - (1.96 * results$s1_dummy)
    results$b1_dummy_uci <- results$b1_dummy + (1.96 * results$s1_dummy)
    results$b2_dummy_lci <- results$b2_dummy - (1.96 * results$s2_dummy)
    results$b2_dummy_uci <- results$b2_dummy + (1.96 * results$s2_dummy)

    results$b1_osca_lci <- results$b1_osca - (1.96 * results$s1_osca)
    results$b1_osca_uci <- results$b1_osca + (1.96 * results$s1_osca)
    results$b2_osca_lci <- results$b2_osca - (1.96 * results$s2_osca)
    results$b2_osca_uci <- results$b2_osca + (1.96 * results$s2_osca)

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
    b1_osca <- results %>% 
        dplyr::group_by(phi) %>%
        dplyr::summarize(binom.test(sum(b1_osca_lci <= v1_mean & b1_osca_uci >= v1_mean), n()) %>% tidy) %>%
        dplyr::select(phi, estimate, conf.low, conf.high) %>%
        dplyr::mutate(genotype="SNP=1", method="OSCA-BF")
    b2_osca <- results %>% 
        dplyr::group_by(phi) %>%
        dplyr::summarize(binom.test(sum(b2_osca_lci <= v2_mean & b2_osca_uci >= v2_mean), n()) %>% tidy) %>%
        dplyr::select(phi, estimate, conf.low, conf.high) %>%
        dplyr::mutate(genotype="SNP=2", method="OSCA-BF")

    coverage1 <- rbind(b1_dummy, b1_osca)
    coverage1 <- merge(coverage1, v1_mean, "phi")
    coverage2 <- rbind(b2_dummy, b2_osca)
    coverage2 <- merge(coverage2, v2_mean, "phi")
    coverage <- rbind(coverage1, coverage2)
    coverage$phi <- as.factor(coverage$phi)
    coverage$genotype <- as.factor(coverage$genotype)

    return(coverage)
}

plot_est <- function(ci){
    p <- ggplot(data=ci, aes(x=x, y=estimate, ymin=conf.low, ymax=conf.high)) +
        geom_point() + 
        geom_abline(intercept = 0, slope = 1, linetype="dashed", color="grey") +
        geom_errorbar(width=.05) +
        theme_classic() +
        xlab("Difference in variance compared with G=0") +
        ylab("Difference in variance compared with G=0 (95% CI)") +
        facet_wrap(method ~ genotype, scales="free") +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
        theme(
            strip.background = element_blank(),
            legend.title.align=0.5,
            text = element_text(size=18)
        )
    return(p)
}

plot_coverage <- function(coverage){
    p <- ggplot(data=coverage, aes(x=x, y=estimate, ymin=conf.low, ymax=conf.high)) +
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
            text = element_text(size=18),
            strip.background = element_blank(),
            legend.title.align=0.5
        )
    return(p)
}

# effect size accuracy
p1 <- plot_est(get_est("sim12_i0/results.csv"))
p2 <- plot_est(get_est("sim12_i1/results.csv"))

# coverage
p3 <- plot_coverage(get_coverage("sim12_i0/results.csv"))
p4 <- plot_coverage(get_coverage("sim12_i1/results.csv"))

# produce multiplot
pdf("sim12.pdf", height=14, width=14)
ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2,  font.label = list(size = 22))
dev.off()