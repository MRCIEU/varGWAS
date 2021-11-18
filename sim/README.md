# Simulations

## Sim1 - Power variance tests to detect change in variance from interaction effect

Simulation

```shell
# normal dist
sbatch runR.sh sim1.R --dist "n"
# mixed normal dist
sbatch runR.sh sim1.R --dist "mn"
# t4 dist
sbatch runR.sh sim1.R --dist "t"
# lognormal dist
sbatch runR.sh sim1.R --dist "l"
```

Plot

```shell
Rscript sim1_plot.R
```

Similar power for Normal data. Higher power for BF/LAD-BF and Levene(mean) for non-Normal. Same power for BF vs LAD-BF.

## Sim2 - T1E of variance tests under no effect with non-normal/normal dist & increasing MAF

Elevated T1E for non-Normal using Levene(mean) and BP. BF/LAD-BF equally well controlled.

## Sim3 - T1E of variance tests under main effect with non-normal/normal dist transformation

Transformations introduce mean-variance effect casusing T1E.

## Sim4 - Runtime performance of OSCA & varGWAS using 10k SNPs and increasing threads

Regression models take 2x longer than non-parametric models. No strong difference between BP vs LAD-BF and Levene vs BF.

## Sim5 - Confounding of the mean and variance effect and adjustment

Adjusting second-stage model for the square of the first-stage model covariates reduces genetic confounding on variance estimate: an example of this - ancestry x SES on T2DM in https://www.thelancet.com/journals/eclinm/article/PIIS2589-5370(21)00240-6/fulltext

## Sim6 - Variance effect estimate and SE

### Deltamethod

sim6b - CIs are correct for var(Y|G==1) but not var(Y|G==2). The latter is too narrow. Although the point estimates are correct for both. Use bootstrap method instead.

### Bootstrap

```shell
# BF-LAD
# perfom reps of sim6
for b in $(seq 0 .5 6); do
    sbatch runR.sh sim6.R -b "$b" -i 1 -n 200
done
# pool reps
echo -n "z " > results.txt; head -n1 results_i1_b0.txt >> results.txt
cat results_i1_b*.txt | grep -v "b1" >> results.txt
Rscript sim6_plot.R
```

CIs have correct coverage using the bootstrap method for var(Y|G==1) and var(Y|G==2) using the dummy method

## Sim7 - FPR of LAD-BF applied to biomarker emperical distribution

Sample with replacement from the emperical biomarker distribution and estimate T1E:

Using the whole distribution - some elevate T1E for highly left-skewed traits

```shell
sbatch runR.sh sim7.R -t "$trait"
```

Dropping observations > +/- 5SD from the mean

```shell
sbatch runR.sh sim7.R -t "$trait" -f
```

Dropping Z >=/- 5SD gave well controlled T1E for all phenotypes

## Sim8 - The OSCA effect estimate

OSCA effect estimate:

- Z-score from P value given normal dist
- Estimate the inverse of the SE given sample size, MAF and Z
- Calculate beta with Z / inverse SE
- Calculate SE with 1/inverse SE
- Update the direction of the beta by estimating the relationship of Y on X

When the outcome is scaled the OSCA effect estimate is ```var(Y)=b/(2/pi)``` a linear relationship between X and var(Y).

Note - scaling the phenotype renders the LAD-BF dummy/x+xsq variance effect estimate invalid.

## Sim9 - Comparing P-val between OSCA - Levene(median), LAD-BF (dummy) and LAD-BF(x+xsq)

OSCA-Levene(median), LAD-BF (dummy) and LAD-BF (x+xsq) give the same P value. Note that dummy/x+xsq is better powered than just having X in the second-stage model even when the effect of X on var(Y) is linear (in which case both models give the correct estimate).

## Sim10 - Variance effect estimate SE comparison between methods

Simulation of linear effect of X on var(Y) i.e. not using an interaction.

Relationship between OSCA effect estimate and true difference in variance is non-linear. SEs are positively correlated for regression and bootstrap models. OSCA SE is inversely correlated with the regression model.

## Sim11 - Second-stage model: X vs X + X^2

Comparison of including X w/wo X^2 in the second-stage model on the estimate for var(Y|G). Having X in the second-stage model allows estimataion when the relationship between X and var(Y) is linear. But an interaction of XU on Y produces a non-linear variance effect of Y conditional on X. Having x+x^2 in the second-stage model or treating X as a dummy variable models the effect correctly.

## Sim12 - Per-genotype effect on var(Y) under interaction

what is var(Y|G==0, G==1, G==2) with both methods? And do the SEs give correct coverage? Also compare with bootstrap

```shell
head -n1 0/sim12_0.csv > results.csv
cat */sim12_*.csv | grep -v b0_dummy >> results.csv
head -n1 0/sim12_0.csv > results.csv
cat */sim12_*.csv | grep -v b0_dummy >> results.csv
```

## Sim13 - Adjusting the variance effect for the interaction

Including U + XU in the second-stage model then the variance effect attenuates