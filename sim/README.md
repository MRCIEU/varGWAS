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

## Sim5 - confounding of the mean and variance effect (not working)

```shell
sbatch runR.sh sim5.R
```

## Sim6 - Variance effect estimate using LAD-BF and OSCA

```shell
# perfom 20 * 20 reps of sim6
for i in {1..20}; do
    sbatch runR.sh sim6.R -b 2 -i "$i" -n 20
done
# pool reps
echo -n "z " > results_b2.txt; head -n1 results_i1_b2.txt >> results_b2.txt
cat results_i*_b2.txt | grep -v "v1" >> results_b2.txt 
```

## Sim7 - false positive rate for subsampled phenotypes

```shell
sbatch runR.sh sim7.R
```