# Simulations

## Sim1 - Power of B-P and Levene's test to detect change in variance from interaction effect

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
Rscript power_plot.R
```

## Sim2 - T1E of B-P and Levene's test under no effect with non-normal/normal dist & increasing MAF

```shell
sbatch runR.sh sim2b.R
```

## Sim3 - T1E of B-P and Levene's test under main effect with non-normal/normal dist transformation

```shell
sbatch runR.sh sim3.R
```

## Sim4 - runtime performance

```shell
sbatch runR.sh sim4.R
```

## Sim5 - confounding of the mean and variance effect

```shell
sbatch runR.sh sim5.R
```

## Sim6 - Linear vs non-linear variance effect estimate

