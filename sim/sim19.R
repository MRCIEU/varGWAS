library("broom")
library("dplyr")
library("ggplot2")
source("funs.R")
set.seed(2)

n_obs <- 1000
n_sim <- 1000
q <- 0.4

# set to have 80% pwr
a <- 0.135
b <- 0.1275
c <- 0.135

results <- data.frame()
for (phi in seq(0, 6, 0.5)){
    for (i in 1:n_sim){
        d <- b * phi
        x <- get_simulated_genotypes(q, n_obs)
        u <- rnorm(n_obs)
        y <- a + b*x + c*u + d*x*u + rnorm(n_obs)
        fit <- lm(y ~ x*u) 
        pi <- fit %>% tidy %>% dplyr::filter(term == "(Intercept)") %>% dplyr::pull(p.value)
        px <- fit %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::pull(p.value)
        pu <- fit %>% tidy %>% dplyr::filter(term == "u") %>% dplyr::pull(p.value)

        yx0 <- var(y[x==0])
        yx1 <- var(y[x==1])
        yx2 <- var(y[x==2])
        
        vu <- 1
        ve <- 1
        A <- c^2*vu+ve # intercept
        B <- 2*c*d*vu # effect of X on var(Y)
        C <- d^2*vu # effect of X^2 on var(Y)
        vy0 <- A+B*0+C*0^2
        vy1 <- A+B*1+C*1^2
        vy2 <- A+B*2+C*2^2

        results <- rbind(results, data.frame(phi, pi, px, pu, gt=0, estimate=yx0, expected=vy0))
        results <- rbind(results, data.frame(phi, pi, px, pu, gt=1, estimate=yx1, expected=vy1))
        results <- rbind(results, data.frame(phi, pi, px, pu, gt=2, estimate=yx2, expected=vy2))
    }
}

# plot results
s <- results %>% 
    dplyr::group_by(phi, gt) %>%
    dplyr::summarize(cbind(expected=expected %>% head(n=1), t.test(estimate, mu=expected %>% head(n=1)) %>% tidy))

s <- s %>% dplyr::filter(gt != 0)
s <- s %>% dplyr::mutate(gt=dplyr::recode(gt, '1'='var(Y|X=1)', '2'='var(Y|X=2)', .default=NA_character_))

pdf("sim19.pdf", width=10)
ggplot(data=s, aes(x=expected, y=estimate, ymin=conf.low, ymax=conf.high, color=phi)) +
    geom_point() + 
    geom_errorbar() +
    theme_classic() + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    geom_abline(xintercept=0, yintercept=1, linetype="dashed", color="grey") +
    xlab("Calculated variance") +
    ylab("Estimated variance") +
    facet_wrap(~gt, scales="free") +
    theme(legend.position = "bottom")
dev.off()