# JLST C++ simulation

## Simulation

```sh
# normal dist
sbatch sim.sh "$n"
# t4 dist
sbatch sim.sh "$t"
# lognormal dist
sbatch sim.sh "$l"
```

## Power

```sh
Rscript power.R --dist "n"
Rscript power.R --dist "l"
Rscript power.R --dist "t"
```