# varGWAS

Implementation of Breusch-Pagan and regression based Brown-Forsythe heteroscedasticity tests for genome-wide association studies

```
T. S. Breusch and A. R. Pagan, “A Simple Test for Heteroscedasticity and Random Coefficient Variation,” Econometrica, vol. 47, no. 5, p. 1287, Sep. 1979, doi: 10.2307/1911963.
M. B. Brown and A. B. Forsythe, “Robust tests for the equality of variances,” J. Am. Stat. Assoc., 1974, doi: 10.1080/01621459.1974.10482955.
```

## Install

Requires UNIX environment

SRC

```shell
git clone git@github.com:MRCIEU/varGWAS.git
cd varGWAS
```

Libraries

```shell
bash lib.sh
```

Build

Tested with GCC v5 & v6. [Libscl](http://www.aronaldg.org/webfiles/libscl/) does not build with GCC >=v7.

```shell
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

Note on HPC systems you may need to explicitly set the compiler path

```shell
# BC4
module load build/gcc-5.5.0
module load tools/cmake/3.20.0
CC=/mnt/storage/software/languages/gcc-5.5.0/bin/gcc \
CXX=/mnt/storage/software/languages/gcc-5.5.0/bin/g++ \
cmake .. -DCMAKE_BUILD_TYPE=Release
```

## Usage

```shell
./varGWAS

Program to perform vGWAS of trait against variants in the BGEN format
Usage:
  varGWAS C++ v1.0.0 [OPTION...]

  -v, --variable_file arg  Path to phenotype file
  -s, --sep arg            File separator
  -c, --covariates arg     List of covariates column names separated by a comma (whitespace and quotes are not permitted).
  -o, --output_file arg    Path to output file
  -b, --bgen_file arg      Path to BGEN file
  -p, --phenotype arg      Column name for phenotype
  -i, --id arg             Column name for genotype identifier
  -r, --robust             Robust method using median value (Brown-Forsythe)
  -h, --help               Print usage
  -t, --threads arg        Number of threads (default: 8)
```

- Unordered categorical variables should be one-hot encoded.
- Do not provide null values in the phenotype file - these should be filtered out.

## Docker

Build image

```shell
docker build -t varGWAS .
```

Perform GWAS

```shell
docker run \
-v /Users/ml18692/projects/varGWAS/test/data:/data \
-e SPDLOG_LEVEL=debug \
-it varGWAS \
-v /data/phenotypes.csv \
-s , \
-c sex,age,PC.1,PC.2,PC.3,PC.4,PC.5,PC.6,PC.7,PC.8,PC.9,PC.10 \
-o /data/output.txt \
-b /data/genotypes.bgen \
-p Y \
-i S \
-t 1
```

## Unit tests

Run unit tests

```shell
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make
```

## Simulation

See [README](./sim/README.md)

## Logging

By default logging level is set to INFO. This can be overidden using environmental variables. See details on
the [spdlog](https://github.com/gabime/spdlog#load-log-levels-from-env-variable-or-from-argv) page.

```shell
export SPDLOG_LEVEL=debug
./varGWAS
```
