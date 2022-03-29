library("dplyr")
library("car")
library("broom")
library("tidyr")
library("GWASTools")
library("jlst")
library("data.table")
library("varGWASR")
library("grid")
library("gtable")
library("ggplot2")
source("funs.R")
set.seed(123)

n_sim <- 1000
n_obs <- 1000

# Sim to evaluate population stratification effects on outcome variance and ability to adjust these effects using LAD-BF vs BF with OLS adjustment

# genotype allele frequency by ancestral groups
AF <- runif(5, min = 0, max = 0.5)

results <- data.frame()
for (A_U1 in c(0, 0.1, 0.2)) {
    for (X_U2 in c(0, 0.01, 0.02)) {
        for (i in 1:n_sim){
            # simulate variables

            # ancestry
            A <- sample(seq(1,5), n_obs, replace = T)

            # genotype
            X <- sapply(A, function(a) {get_simulated_genotypes(AF[a], 1)})
            maf <- mean(X) * .5

            # modifier
            U1 <- rnorm(n_obs)
            U2 <- rnorm(n_obs)

            # outcome
            A <- scale(A)
            X <- scale(X)
            Y <- A*sqrt(.25) + A*U1*sqrt(A_U1) + U1*sqrt(.1) + X*sqrt(.05) + X*U2*sqrt(X_U2) + U2*sqrt(.1) + rnorm(n_obs, sd=sqrt(1-(.25+A_U1+.1+0.05+X_U2+.1)))
            varY <- var(Y)

            data <- data.frame(
                A,
                AU1=A*U1,
                U1,
                X,
                XU2=X*U2,
                U2,
                Y,
                Asq=A^2,
                stringsAsFactors=F
            )

            # OLS regress out confounder
            fit_Y <- lm(Y ~ A, data=data)
            data$Y_A <- resid(fit_Y)

            res <- data.frame(
                rsq.xa=summary(lm(X ~ as.factor(A)))$r.squared,
                rsq.ya=summary(lm(Y ~ A))$r.squared,
                rsq.yau1=summary(lm(Y ~ A*U1))$r.squared,
                rsq.yu1=summary(lm(Y ~ U1))$r.squared,
                rsq.yx=summary(lm(Y ~ X))$r.squared,
                rsq.yxu2=summary(lm(Y ~ X*U2))$r.squared,
                rsq.yu2=summary(lm(Y ~ U2))$r.squared,
                beta_mu0=lm(Y ~ X) %>% tidy %>% dplyr::filter(term=="X") %>% dplyr::pull("estimate"),
                beta_mu1=lm(Y ~ X + A) %>% tidy %>% dplyr::filter(term=="X") %>% dplyr::pull("estimate"),
                p_adj0=vargwas_model(data, "X", "Y", covar1 = NULL, covar2 = NULL) %>% dplyr::pull("phi_p"),
                p_adj1=vargwas_model(data, "X", "Y", covar1 = c("A"), covar2 = NULL) %>% dplyr::pull("phi_p"),
                p_adj2=vargwas_model(data, "X", "Y", covar1 = c("A"), covar2 = c("A")) %>% dplyr::pull("phi_p"),
                p_adj3=vargwas_model(data, "X", "Y", covar1 = c("A"), covar2 = c("Asq")) %>% dplyr::pull("phi_p"),
                p_adj4=vargwas_model(data, "X", "Y", covar1 = c("A"), covar2 = c("A", "Asq")) %>% dplyr::pull("phi_p"),
                p_bf0=leveneTest(data$Y, as.factor(data$X), center=median) %>% tidy %>% dplyr::filter(!is.na(p.value)) %>% dplyr::pull(p.value),
                p_bf1=leveneTest(data$Y_A, as.factor(data$X), center=median) %>% tidy %>% dplyr::filter(!is.na(p.value)) %>% dplyr::pull(p.value),
                A_U1,X_U2,varY,maf
            )

            results <- rbind(results, res)
        }
    }
}

# check R2 is correct
r2 <- results %>%
    dplyr::group_by(A_U1, X_U2) %>%
    dplyr::summarize(
        t.test(rsq.xa) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(rsq.xa="estimate", rsq.xa.low="conf.low", rsq.xa.high="conf.high"),
        t.test(rsq.ya) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(rsq.ya="estimate", rsq.ya.low="conf.low", rsq.ya.high="conf.high"),
        t.test(rsq.yau1) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(rsq.yau1="estimate", rsq.yau1.low="conf.low", rsq.yau1.high="conf.high"),
        t.test(rsq.yu1) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(rsq.yu1="estimate", rsq.yu1.low="conf.low", rsq.yu1.high="conf.high"),
        t.test(rsq.yx) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(rsq.yx="estimate", rsq.yx.low="conf.low", rsq.yx.high="conf.high"),
        t.test(rsq.yxu2) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(rsq.yxu2="estimate", rsq.yxu2.low="conf.low", rsq.yxu2.high="conf.high"),
        t.test(rsq.yu2) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(rsq.yu2="estimate", rsq.yu2.low="conf.low", rsq.yu2.high="conf.high")
    )

# estimate T1E/Power
power <- results %>% 
    dplyr::group_by(A_U1, X_U2) %>%
    dplyr::summarize(
        binom.test(sum(p_adj0 < 0.05), n_sim) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(p_adj0="estimate", p_adj0.low="conf.low", p_adj0.high="conf.high"),
        binom.test(sum(p_adj1 < 0.05), n_sim) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(p_adj1="estimate", p_adj1.low="conf.low", p_adj1.high="conf.high"),
        binom.test(sum(p_adj2 < 0.05), n_sim) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(p_adj2="estimate", p_adj2.low="conf.low", p_adj2.high="conf.high"),
        binom.test(sum(p_adj3 < 0.05), n_sim) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(p_adj3="estimate", p_adj3.low="conf.low", p_adj3.high="conf.high"),
        binom.test(sum(p_adj4 < 0.05), n_sim) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(p_adj4="estimate", p_adj4.low="conf.low", p_adj4.high="conf.high"),
        binom.test(sum(p_bf0 < 0.05), n_sim) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(p_bf0="estimate", p_bf0.low="conf.low", p_bf0.high="conf.high"),
        binom.test(sum(p_bf1 < 0.05), n_sim) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(p_bf1="estimate", p_bf1.low="conf.low", p_bf1.high="conf.high")
    )

long <- rbind(
    power %>% dplyr::select(A_U1, X_U2, p_adj0, p_adj0.low, p_adj0.high) %>% dplyr::rename(estimate="p_adj0", conf.low="p_adj0.low", conf.high="p_adj0.high") %>% dplyr::mutate(estimand = "p_adj0"),
    power %>% dplyr::select(A_U1, X_U2, p_adj1, p_adj1.low, p_adj1.high) %>% dplyr::rename(estimate="p_adj1", conf.low="p_adj1.low", conf.high="p_adj1.high") %>% dplyr::mutate(estimand = "p_adj1"),
    power %>% dplyr::select(A_U1, X_U2, p_adj2, p_adj2.low, p_adj2.high) %>% dplyr::rename(estimate="p_adj2", conf.low="p_adj2.low", conf.high="p_adj2.high") %>% dplyr::mutate(estimand = "p_adj2"),
    power %>% dplyr::select(A_U1, X_U2, p_adj3, p_adj3.low, p_adj3.high) %>% dplyr::rename(estimate="p_adj3", conf.low="p_adj3.low", conf.high="p_adj3.high") %>% dplyr::mutate(estimand = "p_adj3"),
    power %>% dplyr::select(A_U1, X_U2, p_adj4, p_adj4.low, p_adj4.high) %>% dplyr::rename(estimate="p_adj4", conf.low="p_adj4.low", conf.high="p_adj4.high") %>% dplyr::mutate(estimand = "p_adj4"),
    power %>% dplyr::select(A_U1, X_U2, p_bf0, p_bf0.low, p_bf0.high) %>% dplyr::rename(estimate="p_bf0", conf.low="p_bf0.low", conf.high="p_bf0.high") %>% dplyr::mutate(estimand = "p_bf0"),
    power %>% dplyr::select(A_U1, X_U2, p_bf1, p_bf1.low, p_bf1.high) %>% dplyr::rename(estimate="p_bf1", conf.low="p_bf1.low", conf.high="p_bf1.high") %>% dplyr::mutate(estimand = "p_bf1")
)

# plot
long$estimand <- gsub("p_adj0", "LAD-BF_1", long$estimand)
long$estimand <- gsub("p_adj1", "LAD-BF_2", long$estimand)
long$estimand <- gsub("p_adj2", "LAD-BF_3", long$estimand)
long$estimand <- gsub("p_adj3", "LAD-BF_4", long$estimand)
long$estimand <- gsub("p_adj4", "LAD-BF_5", long$estimand)
long$estimand <- gsub("p_bf0", "BF_1", long$estimand)
long$estimand <- gsub("p_bf1", "BF_2", long$estimand)
p <- ggplot(long, aes(x=estimand, y=estimate, ymin=conf.low, ymax=conf.high)) +
    geom_point() +
    geom_errorbar(width=0.3) +
    facet_grid(A_U1~X_U2) +
    theme_classic()+
    labs(y="Proportion P < 0.05 (95% CI)", x="Method") +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey") +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey") +
    ggtitle("Variance explained by genotype interaction effect") +
    theme(
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.spacing = unit(1, "lines"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5, size=11)
    )

# text, size, colour for added text
text = "Variance explained by ancestry interaction effect"
size = 11
col = "black"

# Convert the plot to a grob
gt <- ggplot2::ggplotGrob(p)

# Get the positions of the right strips in the layout: t = top, l = left, ...
strip <-c(subset(gt$layout, grepl("strip-r", gt$layout$name), select = t:r))

# Text grob
text.grob = grid::textGrob(text, rot = -90, gp = gpar(fontsize = size, col = col))

# New column to the right of current strip
# Adjusts its width to text size
width = unit(2, "grobwidth", text.grob) + unit(1, "lines")
gt <- gtable::gtable_add_cols(gt, width, max(strip$r))  

# Add text grob to new column
gt <- gtable::gtable_add_grob(gt, text.grob, 
        t = min(strip$t), l = max(strip$r) + 1, b = max(strip$b))
    

pdf("sim16.pdf")
# Draw it
grid.newpage()
grid.draw(gt)
dev.off()
