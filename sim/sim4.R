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
            Age = rnorm(n_obs),
            Sex = rbinom(n_obs, 2, .5),
            PC1 = rnorm(n_obs),
            PC2 = rnorm(n_obs),
            PC3 = rnorm(n_obs),
            PC4 = rnorm(n_obs),
            PC5 = rnorm(n_obs),
            PC6 = rnorm(n_obs),
            PC7 = rnorm(n_obs),
            PC8 = rnorm(n_obs),
            PC9 = rnorm(n_obs),
            PC10 = rnorm(n_obs),
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
        bf.time <- system.time(system(paste0("varGWAS -v data/phenotypes.csv -s , -o data/gwas-bf.txt -b data/genotypes.bgen -p Y -i S -t ", t)))
        osca_median.time <- system.time(system(paste0("osca --vqtl --bfile data/genotypes --pheno data/phenotypes.txt --out data/osca-median.txt --vqtl-mtd 2 --thread-num ", t)))
        bf.time <- t(data.matrix(bf.time)) %>% as.data.frame
        osca_median.time <- t(data.matrix(osca_median.time)) %>% as.data.frame
        names(bf.time) <- paste0(names(bf.time), ".cpp_bf")
        names(osca_median.time) <- paste0(names(osca_median.time), ".osca_median")

        # store result
        d <- rbind(d, cbind(
            t,
            bf.time,
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
    d %>% dplyr::group_by(t) %>% dplyr::summarize(t.test(elapsed.cpp_bf) %>% tidy %>% select(estimate, conf.low, conf.high) %>% dplyr::mutate(model="LAD Brown-Forsythe", location="Median"))
)
results <- rbind(
    results,
    d %>% dplyr::group_by(t) %>% dplyr::summarize(t.test(elapsed.osca_median) %>% tidy %>% select(estimate, conf.low, conf.high) %>% dplyr::mutate(model="Brown-Forsythe", location="Median"))
)
results$t <- factor(results$t)
results$model <- factor(results$model, levels=c("LAD Brown-Forsythe", "Brown-Forsythe"))

# barchart
pdf("data/sim4a.pdf")
ggplot(data=results, aes(x=model, y=estimate, ymin=conf.low, ymax=conf.high, group=t, color=t, shape=location)) +
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_errorbar(width=.05, position = position_dodge(width = 0.5)) +
    theme_classic() + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    labs(color="Threads", shape="Location") +
    xlab("Method") +
    ylab("Mean runtime (seconds, 95% CI)")
dev.off()