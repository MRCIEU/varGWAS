library("dplyr")
library("broom")
library("tidyr")
library("ggpubr")
library("lmtest")
library("data.table")
source("funs.R")
set.seed(123)

# Requires OSCA and QCTOOL on PATH

n_sim <- 200
n_snps <- 1000
n_obs <- 100000

d <- data.frame()
for (t in c(1,2,4,8)){
    for (i in 1:n_sim){
        # simulate covariates
        data <- data.frame(
            S = paste0("S", seq(1, n_obs)),
            Y = rnorm(n_obs),
            stringsAsFactors=F
        )

        # write out GEN file
        for (j in 1:n_snps){
            write_gen(paste0("data/genotypes.", j, ".gen"), "01", paste0("SNPID_", j), paste0("RSID_", j), j, "A", "G", get_simulated_genotypes(.4, n_obs))
        }

        # combine GEN files
        system("cat data/genotypes.[0-9]*.gen > data/genotypes.gen")

        # write phenotype & sample file
        write.table(file = "data/phenotypes.csv", sep = ",", quote = F, row.names = F, data)
        write.table(file = "data/phenotypes.txt", sep = "\t", quote = F, row.names = F, col.names = F, data[, c("S", "S", "Y")])
        fileConn <- file("data/samples.txt")
        writeLines(c("ID_1 ID_2 missing sex\n0 0 0 D", paste0(data$S, " ", data$S, " ", 0, " ", 1)), fileConn)
        close(fileConn)

        # convert to BGEN file & plink
        system("qctool -g data/genotypes.gen -og data/genotypes.bgen")
        system("bgenix -g data/genotypes.bgen -clobber -index")
        system("qctool -g data/genotypes.gen -s data/samples.txt -og data/genotypes -ofiletype binary_ped")
        system("sed 's/^/S/g' data/genotypes.fam > data/genotypes.fam.sed; mv data/genotypes.fam.sed data/genotypes.fam")

        # run vGWAS
        bp.time <- system.time(system(paste0("varGWAS -v data/phenotypes.csv -s , -o data/gwas-bp.txt -b data/genotypes.bgen -p Y -i S -t ", t)))
        bf.time <- system.time(system(paste0("varGWAS -v data/phenotypes.csv -s , -o data/gwas-bf.txt -b data/genotypes.bgen -p Y -i S -r -t ", t)))
        osca_mean.time <- system.time(system(paste0("osca --vqtl --bfile data/genotypes --pheno data/phenotypes.txt --out data/osca-mean.txt --vqtl-mtd 1 --thread-num ", t)))
        osca_median.time <- system.time(system(paste0("osca --vqtl --bfile data/genotypes --pheno data/phenotypes.txt --out data/osca-median.txt --vqtl-mtd 2 --thread-num ", t)))
        bp.time <- t(data.matrix(bp.time)) %>% as.data.frame
        bf.time <- t(data.matrix(bf.time)) %>% as.data.frame
        osca_mean.time <- t(data.matrix(osca_mean.time)) %>% as.data.frame
        osca_median.time <- t(data.matrix(osca_median.time)) %>% as.data.frame
        names(bp.time) <- paste0(names(bp.time), ".cpp_bp")
        names(bf.time) <- paste0(names(bf.time), ".cpp_bf")
        names(osca_mean.time) <- paste0(names(osca_mean.time), ".osca_mean")
        names(osca_median.time) <- paste0(names(osca_median.time), ".osca_median")

        # store result
        d <- rbind(d, cbind(
            t,
            bp.time,
            bf.time,
            osca_mean.time,
            osca_median.time
        ))
    }
}

# save data
write.table(d, "data/sim4.txt")

# mean and 95% CI of elapsed time
results <- data.frame()
results <- rbind(
    results,
    d %>% dplyr::group_by(t) %>% dplyr::summarize(t.test(elapsed.cpp_bp) %>% tidy %>% select(estimate, conf.low, conf.high) %>% dplyr::mutate(model="Breusch-Pagan", location="Mean"))
)
results <- rbind(
    results,
    d %>% dplyr::group_by(t) %>% dplyr::summarize(t.test(elapsed.cpp_bf) %>% tidy %>% select(estimate, conf.low, conf.high) %>% dplyr::mutate(model="LAD Brown-Forsythe", location="Median"))
)
results <- rbind(
    results,
    d %>% dplyr::group_by(t) %>% dplyr::summarize(t.test(elapsed.osca_mean) %>% tidy %>% select(estimate, conf.low, conf.high) %>% dplyr::mutate(model="Levene", location="Mean"))
)
results <- rbind(
    results,
    d %>% dplyr::group_by(t) %>% dplyr::summarize(t.test(elapsed.osca_median) %>% tidy %>% select(estimate, conf.low, conf.high) %>% dplyr::mutate(model="Brown-Forsythe", location="Median"))
)
results$t <- factor(results$t)
results$model <- factor(results$model, levels=c("Breusch-Pagan", "LAD Brown-Forsythe", "Levene", "Brown-Forsythe"))

# barchart
pdf("data/sim4a.pdf")
ggplot(data=results, aes(x=model, y=estimate, ymin=conf.low, ymax=conf.high, group=t, color=t, shape=location)) +
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_errorbar(width=.05, position = position_dodge(width = 0.5)) +
    theme_classic() + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    labs(color="Location") +
    xlab("Method") +
    ylab("Mean runtime (seconds, 95% CI)")
dev.off()