library("dplyr")
library("broom")
library("tidyr")
library("GWASTools")
library("jlst")
library("data.table")
source("funs.R")
set.seed(123)

n_sim <- 1000
n_obs <- 10000

# Sim to evaluate population stratification effects on outcome variance and ability to adjust these effects using LAD-BF vs BF with OLS adjustment

# genotype allele frequency by ancestral groups
AF <- runif(10, min = 0.05, max = 0.5)

results <- data.frame()
for (i in 1:n_sim){
    # simulate variables

    # ancestry
    A <- sample(seq(1,10), n_obs, replace = T)

    # genotype
    X <- sapply(A, function(a) {get_simulated_genotypes(AF[a], 1)})

    # modifier
    U <- rnorm(n_obs)

    # outcome
    Y <- A*.15 + X*.4 + A*U*.7 + U*.14 + rnorm(n_obs)

    data <- data.frame(
        X,
        Y,
        A,
        U,
        AU=A*U,
        Asq=A^2, Usq=U^2, AUsq=(A*U)^2,
        stringsAsFactors=F
    )

    res <- data.frame(
        rsq.xc=summary(lm(X ~ C))$r.squared,
        rsq.yc=summary(lm(Y ~ C))$r.squared,
        rsq.yx=summary(lm(Y ~ X))$r.squared,
        rsq.ycu=summary(lm(Y ~ C*U))$r.squared,
        beta.yx=lm(Y ~ X) %>% tidy %>% dplyr::filter(term=="X") %>% dplyr::pull("estimate"),
        beta.yxc=lm(Y ~ X + C) %>% tidy %>% dplyr::filter(term=="X") %>% dplyr::pull("estimate"),
        beta.yxcsq=lm(Y ~ X + Csq) %>% tidy %>% dplyr::filter(term=="X") %>% dplyr::pull("estimate"),
        MAF=mean(data$X * .5),
        p_adj1=get_p(data$X, data$Y, covar1 = data %>% dplyr::select("C") %>% as.data.frame, covar2 = NULL),
        p_adj2=get_p(data$X, data$Y, covar1 = data %>% dplyr::select("C") %>% as.data.frame, covar2 = data %>% dplyr::select("C") %>% as.data.frame),
        p_adj3=get_p(data$X, data$Y, covar1 = data %>% dplyr::select("C") %>% as.data.frame, covar2 = data %>% dplyr::select("Csq") %>% as.data.frame),
        p_adj4=get_p(data$X, data$Y, covar1 = data %>% dplyr::select("C") %>% as.data.frame, covar2 = data %>% dplyr::select("C", "Csq") %>% as.data.frame),
        p_adj5=get_p(data$X, data$Y, covar1 = data %>% dplyr::select("C", "Csq") %>% as.data.frame, covar2 = NULL),
        p_adj6=get_p(data$X, data$Y, covar1 = data %>% dplyr::select("C", "Csq") %>% as.data.frame, covar2 = data %>% dplyr::select("C") %>% as.data.frame),
        p_adj7=get_p(data$X, data$Y, covar1 = data %>% dplyr::select("C", "Csq") %>% as.data.frame, covar2 = data %>% dplyr::select("Csq") %>% as.data.frame),
        p_adj8=get_p(data$X, data$Y, covar1 = data %>% dplyr::select("C", "Csq") %>% as.data.frame, covar2 = data %>% dplyr::select("C", "Csq") %>% as.data.frame)

    )

    results <- rbind(results, res)
}

# test power w/wo adjustment
write.csv(results, file="sim5.csv")

# adjusted for C in the first-stage model only 
pdf("data/fs_adj.pdf")
GWASTools::qqPlot(results$p_adj1,main="First-stage model adjusted for C")
dev.off()

# adjusted for C in the first-stage model + C & C^2 in the second-stage
pdf("data/ss_adj.pdf")
GWASTools::qqPlot(results$p_adj4,main="First-stage model adjusted for C & second-stage adjusted for C + C^2")
dev.off()