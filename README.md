# JLSP C++

Implementation of the JLSP method (Staley et al) using C++

```
Staley, J. R., Windmeijer, F., Suderman, M., Smith, G. D., & Tilling, K. (n.d.). A robust mean and variance test with application to epigenome-wide association studies. https://doi.org/10.1101/2020.02.06.926584
```

## Install

Requires UNIX environment

### SRC

```sh
git clone git@ieugit-scmv-d0.epi.bris.ac.uk:ml18692/jlst_cpp.git
cd jlst_cpp
```

### Libraries

```sh
mkdir -p lib
cd lib

# bgen
curl -L http://code.enkre.net/bgen/tarball/release/bgen.tgz > bgen.tar.gz
tar -xvf bgen.tar.gz
mv bgen.tgz bgen
cd bgen
./waf configure
./waf
cd ..

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
curl https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz > eigen-3.3.7.tar.gz
tar -xvf eigen-3.3.7.tar.gz

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

```sh
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

Note on HPC systems you may need to explicitly set the compiler path

```sh
# BC4
module load languages/gcc/9.1.0
module load SQLite/3.13.0-GCC-5.4.0-2.26 

CC=/mnt/storage/software/languages/gcc-9.1/bin/gcc \
CXX=/mnt/storage/software/languages/gcc-9.1/bin/g++ \
../lib/cmake-3.18.6/bin/cmake .. \
-DCMAKE_BUILD_TYPE=Release
```

### Test

TODO - create test bgen file & phenotypes file from a single simulated dataset

Create some test data

```sh
cd test/data
Rscript regression.R
```

Run application tests

```sh
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -G "Unix Makefiles" 
make all
./test/jlst_cpp_test
```

## Usage

```sh
./jlst_cpp
```

Unordered categorical variables should be one-hot encoded.

## Contributing

This project follows the [Google style guide](https://google.github.io/styleguide/cppguide.html)

## Performance

OpenCL/CUDA, OpenMP and MPI?
code optimization and performance tuning, parallelization using both shared memory (OpenMP) and message passing (MPI) paradigms
