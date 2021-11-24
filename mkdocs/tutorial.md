# Tutorial

## Model

The outcome is ```Y = X + U + X*U + E``` where ```X``` is a genotype, ```U``` is a continuous modifier and ```X*U``` is the interaction effect

## Simulate

The script below will simulate the data and requires [qctool](https://www.well.ox.ac.uk/~gav/qctool_v2/) on the PATH.

```shell
Rscript test/data/example.R
```

Alternatively the data are provided in ```test/data```.

## GWAS

Test for the effect of the SNP on the variance of the outcome

```shell
./varGWAS \
-v test/data/example.csv \
-s , \
-o test/data/example.txt \
-b test/data/example.bgen \
-p Y \
-i S
```

## Output

The effect of the SNP on outcome variance is non-linear so the genotype is treated as a dummy variable in the second-stage regression. This means there are two effects of the SNP-var(Y) relationship for each level of the genotype.

| chr | pos | rsid   | oa | ea | n     | eaf     | beta         | se        | t           | p        | theta       | phi_x1   | se_x1     | phi_x2  | se_x2    | phi_f   | phi_p        |
|-----|-----|--------|----|----|-------|---------|--------------|-----------|-------------|----------|-------------|----------|-----------|---------|----------|---------|--------------|
| 01  | 1   | RSID_1 | G  | A  | 10000 | 0.39485 | -0.000127464 | 0.0144545 | -0.00881832 | 0.992964 | -0.00143247 | 0.489362 | 0.0267757 | 1.85565 | 0.095883 | 667.129 | 1.09461e-272 |

- ```chr```, ```pos```, ```rsid```, ```oa``` (non-effect allele) and ```ea``` (effect allele) describe the variant
- ```n``` and ```eaf``` are the total sample size and effect allele frequency included in the model
- ```beta```, ```se```,  ```t``` and ```p``` describe the effect of the SNP on the mean of the outcome
- ```theta``` is the effect of the SNP on the median of the outcome
- ```phi_x1``` and ```phi_x2``` is the average change in variance from ```SNP=0``` to ```SNP=1``` and ```SNP=2```. ```se_x1``` and ```se_x2``` are the standard errors of these statistics.
- ```phi_f``` and ```phi_p``` are the F-statistic and P-value for the effect of the SNP on outcome variance

The trait was standardised (see ```test/data/example.R```) so the units are ```sigma^2```, SNP=1 was associated with an increase of 0.489 and 1.856 for SNP=2.  