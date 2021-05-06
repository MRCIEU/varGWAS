# Simulations

## Power

Simulation

```shell
# normal dist
sbatch runR.sh sim1.R --dist "n"
# t4 dist
sbatch runR.sh sim1.R --dist "t"
# lognormal dist
sbatch runR.sh sim1.R --dist "l"
```

Plot

```shell
Rscript power.R --dist "n"
Rscript power.R --dist "l"
Rscript power.R --dist "t"
```

## Confounding

```shell
sbatch runR.sh sim2.R
```