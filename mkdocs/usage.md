# Usage

Requires UNIX environment. Use [Docker](install.md#Docker) for Windows.

```shell
./varGWAS

Program to perform GWAS of trait variability against variants in the BGEN format
Usage:
  varGWAS v1.2.1 [OPTION...]

  -v, --variable_file arg  Path to phenotype file
  -s, --sep arg            File separator
  -c, --covariates arg     List of covariates column names separated by a comma (whitespace and quotes are not permitted).
  -o, --output_file arg    Path to output file
  -b, --bgen_file arg      Path to BGEN file
  -p, --phenotype arg      Column name for phenotype
  -i, --id arg             Column name for genotype identifier
  -m, --maf arg            Filter out variants with a MAF below this threshold
  -h, --help               Print usage
  -t, --threads arg        Number of threads
```

## Phenotypes

- Do not provide null values in the phenotype file - these should be filtered out.
- Unordered categorical variables should be one-hot encoded (dummy variables).
- Include the square of continuous/ordinal phenotypes to adjust the variance effect.
- The variance effect size is a unitless measure; standardise the outcome beforehand by dividing the trait by its SD.

## Output

See description of GWAS summary stats [here](tutorial.md#Output) 

# Logging

By default logging level is set to INFO. This can be overidden using environmental variables. See details on
the [spdlog](https://github.com/gabime/spdlog#load-log-levels-from-env-variable-or-from-argv) page.

```shell
export SPDLOG_LEVEL=debug
./varGWAS
```

# Unit tests

Run unit tests

```shell
# build debug release
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make
# run tests
./bin/varGWAS_test
```

# Simulations

See [README](https://github.com/MRCIEU/varGWAS/blob/master/sim/README.md) for simulations of test power, type 1 error, accuracy and coverage etc.

# Issues

Report issues [here](https://github.com/MRCIEU/varGWAS/issues)