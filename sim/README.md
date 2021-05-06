# Simulations

## Power

Simulation

```sh
# normal dist
sbatch runR.sh sim1.R --dist "n"
# t4 dist
sbatch runR.sh sim1.R --dist "t"
# lognormal dist
sbatch runR.sh sim1.R --dist "l"
```

Plot

```sh
Rscript power.R --dist "n"
Rscript power.R --dist "l"
Rscript power.R --dist "t"
```

## Confounding
