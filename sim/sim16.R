library("dplyr")
library("broom")
library("tidyr")
library("GWASTools")
library("jlst")
library("data.table")
library("varGWASR")
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
    A <- scale(A)
    X <- scale(X)
    Y <- A*sqrt(.2) + X*sqrt(.05) + A*U*sqrt(.1) + U*sqrt(.2) + rnorm(n_obs, sd=sqrt(1-(.2+.05+.1+.2)))
    varY <- var(Y)

    # OLS regress out confounder
    fit_Y <- lm(Y ~ A, data=data)
    fit_X <- lm(X ~ A, data=data)

    data <- data.frame(
        X,
        Y,
        Y_A=resid(fit_Y),
        X_A=resid(fit_X),
        A,
        U,
        AU=A*U,
        Asq=A^2, Usq=U^2, AUsq=(A*U)^2,
        stringsAsFactors=F
    )

    res <- data.frame(
        varY,
        rsq.xa=summary(lm(X ~ as.factor(A)))$r.squared,
        rsq.ya=summary(lm(Y ~ A))$r.squared,
        rsq.yx=summary(lm(Y ~ X))$r.squared,
        rsq.yau=summary(lm(Y ~ A*U))$r.squared,
        beta_mu0=lm(Y ~ X) %>% tidy %>% dplyr::filter(term=="X") %>% dplyr::pull("estimate"),
        beta_mu1=lm(Y ~ X + A) %>% tidy %>% dplyr::filter(term=="X") %>% dplyr::pull("estimate"),
        p_adj0=vargwas_model(data, "X", "Y", covar1 = NULL, covar2 = NULL) %>% dplyr::pull("phi_p"),
        p_adj1=vargwas_model(data, "X", "Y", covar1 = c("A"), covar2 = NULL) %>% dplyr::pull("phi_p"),
        p_adj2=vargwas_model(data, "X", "Y", covar1 = c("A"), covar2 = c("A")) %>% dplyr::pull("phi_p"),
        p_adj3=vargwas_model(data, "X", "Y", covar1 = c("A"), covar2 = c("Asq")) %>% dplyr::pull("phi_p"),
        p_adj4=vargwas_model(data, "X", "Y", covar1 = c("A"), covar2 = c("A", "Asq")) %>% dplyr::pull("phi_p"),
        p_bf0=leveneTest(data$Y, as.factor(data$X), center=median) %>% tidy %>% dplyr::filter(!is.na(p.value)) %>% dplyr::pull(p.value),
        p_bf1=leveneTest(data$Y_A, as.factor(data$X), center=median) %>% tidy %>% dplyr::filter(!is.na(p.value)) %>% dplyr::pull(p.value)
    )

    results <- rbind(results, res)
}

# estimate T1E
results %>% dplyr::summarize(
    binom.test(sum(p_adj0 < 0.05), n_sim) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(p_adj0="estimate", p_adj0.low="conf.low", p_adj0.high="conf.high"),
    binom.test(sum(p_adj1 < 0.05), n_sim) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(p_adj1="estimate", p_adj1.low="conf.low", p_adj1.high="conf.high"),
    binom.test(sum(p_adj2 < 0.05), n_sim) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(p_adj2="estimate", p_adj2.low="conf.low", p_adj2.high="conf.high"),
    binom.test(sum(p_adj3 < 0.05), n_sim) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(p_adj3="estimate", p_adj3.low="conf.low", p_adj3.high="conf.high"),
    binom.test(sum(p_adj4 < 0.05), n_sim) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(p_adj4="estimate", p_adj4.low="conf.low", p_adj4.high="conf.high"),
    binom.test(sum(p_bf0 < 0.05), n_sim) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(p_bf0="estimate", p_bf0.low="conf.low", p_bf0.high="conf.high"),
    binom.test(sum(p_bf1 < 0.05), n_sim) %>% tidy %>% dplyr::select(estimate, conf.low, conf.high) %>% dplyr::rename(p_bf1="estimate", p_bf1.low="conf.low", p_bf1.high="conf.high")
)