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

## Sim4 - Runtime performance of OSCA & varGWAS

```shell
sbatch runR.sh sim4.R
```

## Sim5 - Confounding of the mean and variance effect (TODO)

```shell
sbatch runR.sh sim5.R
```

## Sim6 - Variance effect bootstrap

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

## Sim7 - false positive rate for subsampled phenotypes

```shell
sbatch runR.sh sim7.R -t "$trait" -f
```

## Sim8 - variance effect estimate

## Sim9 - P value comparison

OSCA - BF, LAD-BF (dummy) and LAD-BF (xsq) give the same P value
Check the OSCA method to derive the BETA and SE

```shell
sbatch runR.sh sim9.R
```

## Sim10 - SE comparison

```shell
sbatch runR.sh sim10.R
```

## Sim11 - X vs x + sq

