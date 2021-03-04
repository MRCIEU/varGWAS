# JLSP C++

Implementation of the JLSP method (Staley et al) using C++

```
Staley, J. R., Windmeijer, F., Suderman, M., Smith, G. D., & Tilling, K. (n.d.). A robust mean and variance test with application to epigenome-wide association studies. https://doi.org/10.1101/2020.02.06.926584
```

## Install

Requires UNIX environment

### SRC

```shell
git clone git@ieugit-scmv-d0.epi.bris.ac.uk:ml18692/jlst_cpp.git
cd jlst_cpp
```

### Libraries

```shell
bash lib.sh
```

### Build

Tested with GCC v5 & v6 and with Apple Clang v12. Quantile regression library does not build with GCC >=v7.

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
PATH="$PATH":/mnt/storage/home/ml18692/projects/jlst_cpp/lib/cmake-3.20.0-rc2-linux-x86_64/bin
CMAKE_ROOT=/mnt/storage/home/ml18692/projects/jlst_cpp/lib/cmake-3.20.0-rc2-linux-x86_64
CC=/mnt/storage/software/languages/gcc-5.5.0/bin/gcc \
CXX=/mnt/storage/software/languages/gcc-5.5.0/bin/g++ \
cmake .. \
-DCMAKE_BUILD_TYPE=Debug
```

### Docker

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

### Test

Run unit tests

```shell
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make
```

Run end-to-end test

```shell
cd test
Rscript sim.R
```

Inspect output [results.csv](./test/data/results.csv)

## Usage

```shell
./bin/jlst_cpp

Program to perform vGWAS of trait against variants in the BGEN format
Usage:
  JLST C++ v0.0.1 [OPTION...]

  -v, --variable_file arg  Path to phenotype file
  -s, --sep arg            File separator
  -c, --covariates arg     List of covariates column names separated by a comma (whitespace and quotes are not permitted).
  -o, --output_file arg    Path to output file
  -b, --bgen_file arg      Path to BGEN file
  -p, --phenotype arg      Column name for phenotype
  -i, --id arg             Column name for genotype identifier
  -h, --help               Print usage
  -t, --threads arg        Number of threads (default: 8)
```

- Unordered categorical variables should be one-hot encoded.
- Do not provide null values in the phenotype file - these should be filtered out.

## Logging

By default logging level is set to INFO. This can be overidden using environmental variables. See details on
the [spdlog](https://github.com/gabime/spdlog#load-log-levels-from-env-variable-or-from-argv) page.

```shell
export SPDLOG_LEVEL=debug
./bin/jlst_cpp
```

## Contributing

This project follows the [Google style guide](https://google.github.io/styleguide/cppguide.html)

## Performance

OpenCL/CUDA, OpenMP and MPI? code optimization and performance tuning, parallelization using both shared memory (OpenMP)
and message passing (MPI) paradigms