# JLSP C++

Implementation of the JLSP method (Staley et al) using C++

```
Staley, J. R., Windmeijer, F., Suderman, M., Smith, G. D., & Tilling, K. (n.d.). A robust mean and variance test with application to epigenome-wide association studies. https://doi.org/10.1101/2020.02.06.926584
```

## Install

### SRC

```sh
git clone git@ieugit-scmv-d0.epi.bris.ac.uk:ml18692/jlst_cpp.git
cd jlst_cpp
```

### Libraries

```sh
mkdir -p includes
cd includes

# bgen
wget http://bitbucket.org/gavinband/bgen/get/master.tar.gz
tar -xvf master.tar.gz

# boost
wget https://dl.bintray.com/boostorg/release/1.72.0/source/boost_1_72_0.tar.gz
tar -xvf boost_1_72_0.tar.gz
cd boost_1_72_0
mkdir build
./bootstrap.sh
./b2 --with-test --prefix=$PWD/build install

# google test
wget https://github.com/google/googletest/archive/release-1.10.0.tar.gz
tar -xvf release-1.10.0.tar.gz
mkdir tests/lib
mv googletest-release-1.10.0 tests/lib

# quantile regression
wget --no-parent -r http://www.aronaldg.org/webfiles/compecon/src/libscl_float/
wget --no-parent -r http://www.aronaldg.org/webfiles/compecon/src/libscl_demo/
cp libscl_float/src/scltypes.tpl libscl_float/src/scltypes.h

# linear regression
wget https://www.mlpack.org/files/mlpack-3.2.2.tar.gz
tar -xvf mlpack-3.2.2.tar.gz

cd ..
```

### Build

```sh
mkdir -p build
cd build
cmake ..
make
```

## Usage

```sh
./jlst_cpp ../includes/bgen/example/example.v11.bgen default
```
