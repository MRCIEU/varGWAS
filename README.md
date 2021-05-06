# JLSP C++

Implementation of the Breusch-Pagan heteroscedasticity test for genome-wide association studies

```
T. S. Breusch and A. R. Pagan, “A Simple Test for Heteroscedasticity and Random Coefficient Variation,” Econometrica, vol. 47, no. 5, p. 1287, Sep. 1979, doi: 10.2307/1911963.
```

## Install

Requires UNIX environment

SRC

```shell
git clone git@ieugit-scmv-d0.epi.bris.ac.uk:ml18692/jlst_cpp.git
cd jlst_cpp
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
./bin/jlst_cpp

Program to perform vGWAS of trait against variants in the BGEN format
Usage:
  JLST C++ v1.0.0 [OPTION...]

  -v, --variable_file arg  Path to phenotype file
  -s, --sep arg            File separator
  -c, --covariates arg     List of covariates column names separated by a comma (whitespace and quotes are not permitted).
  -o, --output_file arg    Path to output file
  -b, --bgen_file arg      Path to BGEN file
  -p, --phenotype arg      Column name for phenotype
  -i, --id arg             Column name for genotype identifier
  -h, --help               Print usage
  -t, --threads arg        Number of threads
```

- Unordered categorical variables should be one-hot encoded.
- Do not provide null values in the phenotype file - these should be filtered out.

## Docker

Build image

```shell
docker build -t jlst_cpp .
```

Perform vGWAS

```shell
docker run \
-v /Users/ml18692/projects/jlst_cpp/test/data:/data \
-e SPDLOG_LEVEL=debug \
-it jlst_cpp \
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
./bin/jlst_cpp
```

## Contributing

This project follows the [Google style guide](https://google.github.io/styleguide/cppguide.html)