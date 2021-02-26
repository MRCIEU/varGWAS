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
mkdir -p lib
cd lib

# boost
curl -L https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.bz2 > boost_1_75_0.tar.bz2
tar --bzip2 -xf boost_1_75_0.tar.bz2

# bgen
curl -L http://code.enkre.net/bgen/tarball/release/bgen.tgz > bgen.tar.gz
tar -xvf bgen.tar.gz
mv bgen.tgz bgen
cd bgen
./waf configure
./waf
cd ..

# zstd
curl -L https://github.com/facebook/zstd/archive/v1.1.0.tar.gz > v1.1.0.tar.gz
tar -xvf v1.1.0.tar.gz
cd zstd-1.1.0/build/cmake
mkdir build
cd build
cmake ..
make

# quantile regression
curl -L http://www.aronaldg.org/webfiles/libscl/libscl.tar > libscl.tar
tar -xvf libscl.tar
cd libscl/gpp
make
cd ../..

# csv reader 
curl -L https://github.com/ben-strasser/fast-cpp-csv-parser/archive/master.zip > fast-cpp-csv-parser.zip
unzip fast-cpp-csv-parser.zip
mv fast-cpp-csv-parser-master fast-cpp-csv-parser

# eigen
curl https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz > eigen-3.3.9.tar.gz
tar -xvf eigen-3.3.9.tar.gz

# cxxopts
curl -L https://github.com/jarro2783/cxxopts/archive/v2.2.0.tar.gz > cxxopts-2.2.0.tar.gz
tar -xvf cxxopts-2.2.0.tar.gz

# spdlog
git clone https://github.com/gabime/spdlog.git
cd spdlog && mkdir build && cd build
cmake .. && make -j

# google test
curl -L https://github.com/google/googletest/archive/release-1.10.0.tar.gz > release-1.10.0.tar.gz
tar -xvf release-1.10.0.tar.gz

# ThreadPool
curl -L https://github.com/progschj/ThreadPool/archive/master.zip > ThreadPool.zip
unzip ThreadPool.zip
mv ThreadPool-master ThreadPool
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

Unordered categorical variables should be one-hot encoded.
Do not provide null values in the phenotype file - these should be filtered out.

## Contributing

This project follows the [Google style guide](https://google.github.io/styleguide/cppguide.html)

## Performance

OpenCL/CUDA, OpenMP and MPI?
code optimization and performance tuning, parallelization using both shared memory (OpenMP) and message passing (MPI) paradigms