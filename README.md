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
CC=/mnt/storage/software/languages/gcc-5.5.0/bin/gcc \
CXX=/mnt/storage/software/languages/gcc-5.5.0/bin/g++ \
../lib/cmake-3.18.6/bin/cmake .. \
-DCMAKE_BUILD_TYPE=Release
```

### Docker

```shell
docker build -t jlst_cpp .
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