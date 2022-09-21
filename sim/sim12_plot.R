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
    b1_drm_ci <- results %>%
        dplyr::group_by(phi) %>%
        dplyr::summarize(t.test(b1_drm) %>% tidy) %>%
        dplyr::select(phi, estimate, conf.low, conf.high) %>%
        dplyr::mutate(genotype="SNP=1", method="DRM")
    b2_drm_ci <- results %>%
        dplyr::group_by(phi) %>%
        dplyr::summarize(t.test(b2_drm) %>% tidy) %>%
        dplyr::select(phi, estimate, conf.low, conf.high) %>%
        dplyr::mutate(genotype="SNP=2", method="DRM")
    b1_quail_ci <- results %>%
        dplyr::group_by(phi) %>%
        dplyr::summarize(t.test(b1_quail) %>% tidy) %>%
        dplyr::select(phi, estimate, conf.low, conf.high) %>%
        dplyr::mutate(genotype="SNP=1", method="QUAIL")
    b2_quail_ci <- results %>%
        dplyr::group_by(phi) %>%
        dplyr::summarize(t.test(b2_quail) %>% tidy) %>%
        dplyr::select(phi, estimate, conf.low, conf.high) %>%
        dplyr::mutate(genotype="SNP=2", method="QUAIL")
    v1_mean <- results %>%
        dplyr::group_by(phi) %>%
        dplyr::summarize(x=mean(v1-v0))
    v2_mean <- results %>%
        dplyr::group_by(phi) %>%
        dplyr::summarize(x=mean(v2-v0))
    
    # rescale estimates
    b1_osca_ci$estimate <- b1_osca_ci$estimate / (2/pi)
    b2_osca_ci$estimate <- b2_osca_ci$estimate / (2/pi)
    b1_osca_ci$conf.low <- b1_osca_ci$conf.low / (2/pi)
    b1_osca_ci$conf.high <- b1_osca_ci$conf.high / (2/pi)
    b2_osca_ci$conf.low <- b2_osca_ci$conf.low / (2/pi)
    b2_osca_ci$conf.high <- b2_osca_ci$conf.high / (2/pi)

    b1_drm_ci$estimate <- b1_drm_ci$estimate * sqrt(2*pi)
    b2_drm_ci$estimate <- b2_drm_ci$estimate * sqrt(2*pi)
    b1_drm_ci$conf.low <- b1_drm_ci$conf.low * sqrt(2*pi)
    b1_drm_ci$conf.high <- b1_drm_ci$conf.high * sqrt(2*pi)
    b2_drm_ci$conf.low <- b2_drm_ci$conf.low * sqrt(2*pi)
    b2_drm_ci$conf.high <- b2_drm_ci$conf.high * sqrt(2*pi)

    b1_quail_ci$estimate <- b1_quail_ci$estimate / sqrt(2/pi)
    b2_quail_ci$estimate <- b2_quail_ci$estimate / sqrt(2/pi)
    b1_quail_ci$conf.low <- b1_quail_ci$conf.low / sqrt(2/pi)
    b1_quail_ci$conf.high <- b1_quail_ci$conf.high / sqrt(2/pi)
    b2_quail_ci$conf.low <- b2_quail_ci$conf.low / sqrt(2/pi)
    b2_quail_ci$conf.high <- b2_quail_ci$conf.high / sqrt(2/pi)

    v1 <- rbind(b1_dummy_ci, b1_osca_ci, b1_drm_ci, b1_quail_ci)
    v1 <- merge(v1, v1_mean, "phi")
    v2 <- rbind(b2_dummy_ci, b2_osca_ci, b2_drm_ci, b2_quail_ci)
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
    names(results)[23] <- "v1_mean"
    names(results)[24] <- "v2_mean"
    results$b1_dummy_lci <- results$b1_dummy - (1.96 * results$s1_dummy)
    results$b1_dummy_uci <- results$b1_dummy + (1.96 * results$s1_dummy)
    results$b2_dummy_lci <- results$b2_dummy - (1.96 * results$s2_dummy)
    results$b2_dummy_uci <- results$b2_dummy + (1.96 * results$s2_dummy)

    results$b1_osca_lci <- results$b1_osca - (1.96 * results$s1_osca)
    results$b1_osca_uci <- results$b1_osca + (1.96 * results$s1_osca)
    results$b2_osca_lci <- results$b2_osca - (1.96 * results$s2_osca)
    results$b2_osca_uci <- results$b2_osca + (1.96 * results$s2_osca)

    results$b1_drm_lci <- results$b1_drm - (1.96 * results$s1_drm)
    results$b1_drm_uci <- results$b1_drm + (1.96 * results$s1_drm)
    results$b2_drm_lci <- results$b2_drm - (1.96 * results$s2_drm)
    results$b2_drm_uci <- results$b2_drm + (1.96 * results$s2_drm)

    results$b1_quail_lci <- results$b1_quail - (1.96 * results$s1_quail)
    results$b1_quail_uci <- results$b1_quail + (1.96 * results$s1_quail)
    results$b2_quail_lci <- results$b2_quail - (1.96 * results$s2_quail)
    results$b2_quail_uci <- results$b2_quail + (1.96 * results$s2_quail)

    # rescale estimates
    results$b1_osca <- results$b1_osca / (2/pi)
    results$b2_osca <- results$b2_osca / (2/pi)
    results$b1_osca_lci <- results$b1_osca_lci / (2/pi)
    results$b1_osca_uci <- results$b1_osca_uci / (2/pi)
    results$b2_osca_lci <- results$b2_osca_lci / (2/pi)
    results$b2_osca_uci <- results$b2_osca_uci / (2/pi)

    results$b1_drm <- results$b1_drm * sqrt(2*pi)
    results$b2_drm <- results$b2_drm * sqrt(2*pi)
    results$b1_drm_lci <- results$b1_drm_lci * sqrt(2*pi)
    results$b1_drm_uci <- results$b1_drm_uci * sqrt(2*pi)
    results$b2_drm_lci <- results$b2_drm_lci * sqrt(2*pi)
    results$b2_drm_uci <- results$b2_drm_uci * sqrt(2*pi)

    results$b1_quail <- results$b1_quail / sqrt(2/pi)
    results$b2_quail <- results$b2_quail / sqrt(2/pi)
    results$b1_quail_lci <- results$b1_quail_lci / sqrt(2/pi)
    results$b1_quail_uci <- results$b1_quail_uci / sqrt(2/pi)
    results$b2_quail_lci <- results$b2_quail_lci / sqrt(2/pi)
    results$b2_quail_uci <- results$b2_quail_uci / sqrt(2/pi)

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
    b1_drm <- results %>% 
        dplyr::group_by(phi) %>%
        dplyr::summarize(binom.test(sum(b1_drm_lci <= v1_mean & b1_drm_uci >= v1_mean), n()) %>% tidy) %>%
        dplyr::select(phi, estimate, conf.low, conf.high) %>%
        dplyr::mutate(genotype="SNP=1", method="DRM")
    b2_drm <- results %>% 
        dplyr::group_by(phi) %>%
        dplyr::summarize(binom.test(sum(b2_drm_lci <= v2_mean & b2_drm_uci >= v2_mean), n()) %>% tidy) %>%
        dplyr::select(phi, estimate, conf.low, conf.high) %>%
        dplyr::mutate(genotype="SNP=2", method="DRM")
    b1_quail <- results %>% 
        dplyr::group_by(phi) %>%
        dplyr::summarize(binom.test(sum(b1_quail_lci <= v1_mean & b1_quail_uci >= v1_mean), n()) %>% tidy) %>%
        dplyr::select(phi, estimate, conf.low, conf.high) %>%
        dplyr::mutate(genotype="SNP=1", method="QUAIL")
    b2_quail <- results %>% 
        dplyr::group_by(phi) %>%
        dplyr::summarize(binom.test(sum(b2_quail_lci <= v2_mean & b2_quail_uci >= v2_mean), n()) %>% tidy) %>%
        dplyr::select(phi, estimate, conf.low, conf.high) %>%
        dplyr::mutate(genotype="SNP=2", method="QUAIL")

    coverage1 <- rbind(b1_dummy, b1_osca, b1_drm, b1_quail)
    coverage1 <- merge(coverage1, v1_mean, "phi")
    coverage2 <- rbind(b2_dummy, b2_osca, b2_drm, b2_quail)
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
        facet_wrap(method ~ genotype, scales="free",  ncol = 2) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
        theme(
            strip.background = element_blank(),
            legend.title.align=0.5
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
            strip.background = element_blank(),
            legend.title.align=0.5
        )
    return(p)
}

# effect size accuracy
p1 <- plot_est(get_est("sim12_i0_new/results.csv"))
p2 <- plot_est(get_est("sim12_i1_new/results.csv"))

# coverage
p3 <- plot_coverage(get_coverage("sim12_i0_new/results.csv"))
p4 <- plot_coverage(get_coverage("sim12_i1_new/results.csv"))

# print plots
pdf("sim12_p1.pdf")
print(p1)
dev.off()
pdf("sim12_p2.pdf")
print(p2)
dev.off()
pdf("sim12_p3.pdf")
print(p3)
dev.off()
pdf("sim12_p4.pdf")
print(p4)
dev.off()