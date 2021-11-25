library("dplyr")
library("broom")
library("tidyr")
library("cqrReg")
library("GWASTools")
library("jlst")
source("funs.R")
library("data.table")
set.seed(123)

n_sim <- 1000
n_obs <- 10000

# get P value using LAD-BF
get_p <- function(x, y, covar1=NULL, covar2=NULL){
    if (!is.null(covar1)){
        X <- as.matrix(cbind(x, covar1))
    } else {
        X <- as.matrix(data.frame(x))
    }
    # first-stage fit
    fit <- qrfit(X=X, y=y, tau=.5, method="mm")
    b <- rbind(fit$b, fit$beta)
    # predicted
    X <- cbind(rep(1, nrow(X)), X)
    fitted <- X %*% b
    # residual
    d <- y - fitted
    # abs residual
    d <- abs(as.vector(d))
    # second-stage model
    x <- as.factor(x)
    if (!is.null(covar2)){
        X <- data.frame(x)
        X <- cbind(X, covar2)
        fit_null <- lm(d ~ ., data=covar2)
    } else {
        X <- data.frame(x)
        fit_null <- lm(d ~ 1)
    }
    fit2 <- lm(d ~ ., data=X)
    p <- anova(fit_null, fit2) %>% tidy %>% dplyr::pull(p.value) %>% dplyr::nth(2)
    return(p)
}

results <- data.frame()
for (i in 1:n_sim){
    # simulate variables

    # confounder
    C <- rnorm(n_obs)
    Csq <- C^2
    # modifier
    U <- rnorm(n_obs)
    Usq <- U^2
    # interaction
    CU <- C*U
    CUsq <- CU^2
    # genotype with MAF of 0.4
    # confounder explains 20% of X
    X <- sapply(C, function(c) {p <- 1 / (1 + exp(-c)); get_simulated_genotypes(p * .8, 1)})
    # confounder explains 5% of Y
    # X explains 7.5% of Y
    # CU explains 40% of Y
    # outcome
    Y <- C*.15 + X*.4 + C*U*.7 + rnorm(n_obs)

    data <- data.frame(
        X,
        Y,
        C,
        U,
        CU,
        Csq, Usq, CUsq,
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