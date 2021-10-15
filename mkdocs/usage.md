# Usage

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

# Docker

Perform GWAS

```shell
docker run \
-v /Users/ml18692/projects/varGWAS/test/data:/data \
-e SPDLOG_LEVEL=debug \
-it vargwas \
-v /data/phenotypes.csv \
-s , \
-c sex,age,PC.1,PC.2,PC.3,PC.4,PC.5,PC.6,PC.7,PC.8,PC.9,PC.10 \
-o /data/output.txt \
-b /data/genotypes.bgen \
-p Y \
-i S \
-t 1
```

# Unit tests

Run unit tests

```shell
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make
```

# Simulation

See [README](https://github.com/MRCIEU/varGWAS/blob/master/sim/README.md)

# Logging

By default logging level is set to INFO. This can be overidden using environmental variables. See details on
the [spdlog](https://github.com/gabime/spdlog#load-log-levels-from-env-variable-or-from-argv) page.

```shell
export SPDLOG_LEVEL=debug
./varGWAS
```
