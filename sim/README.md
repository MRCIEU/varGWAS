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
Rscript power_plot.R
```

## T1E

```shell
sbatch runR.sh sim2.R
```

Plot

```shell
Rscript t1e_plot.R
```