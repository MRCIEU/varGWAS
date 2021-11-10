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
Rscript sim1_plot.R
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
# perfom reps of sim6
for b in $(seq 0 .5 6); do
    sbatch runR.sh sim6.R -b "$b" -i 1 -n 200
done
# pool reps
echo -n "z " > results.txt; head -n1 results_i1_b0.txt >> results.txt
cat results_i1_b*.txt | grep -v "b1" >> results.txt
sbatch runR.sh sim6b.R
```

## Sim7 - false positive rate for subsampled phenotypes

```shell
sbatch runR.sh sim7.R
```