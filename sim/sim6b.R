library("data.table")
library("dplyr")
library("broom")
library("ggplot2")
set.seed(13)

# load results
d <- fread("results.txt")

# calculate expresssion for true variance difference for each value of beta
t <- d %>% 
    dplyr::group_by(b) %>%
    dplyr::summarize(t1_t=mean(t1), t2_t=mean(t2))
t$t1_t <- 1^2*t$b^2
t$t2_t <- 2^2*t$b^2

# append true variance difference to bs results
d <- merge(d, t, "b")

# calculate coverage of bootstrap replicates & check CI includes 95%
r1 <- d %>% 
    dplyr::group_by(b) %>%
    dplyr::summarize(binom.test(sum(lci1 <= t1_t & uci1 >= t1_t), n()) %>% tidy)
r2 <- d %>% 
    dplyr::group_by(b) %>%
    dplyr::summarize(binom.test(sum(lci2 <= t2_t & uci2 >= t2_t), n()) %>% tidy)
r1 %>% dplyr::filter(conf.low > .95 | conf.high < .95)
r2 %>% dplyr::filter(conf.low > .95 | conf.high < .95)

# plot expected vs measured variance differences
ggplot(data=d, aes(x=t2_t, y=b2, ymin=lci2, ymax=uci2)) +
    geom_point() + 
    geom_errorbar(width=.05) +
    theme_classic() + 
    xlab("True difference in variance between SNP=0 and SNP=2") +
    ylab("Estimated difference in variance") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))