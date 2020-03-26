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
curl https://bitbucket.org/gavinband/bgen/get/44fcabbc5c38.zip > gavinband-bgen-44fcabbc5c38.zip
unzip gavinband-bgen-44fcabbc5c38.zip
cd gavinband-bgen-44fcabbc5c38
./waf configure
./waf
cd ..

# quantile regression
wget http://www.aronaldg.org/webfiles/libscl/libscl.tar
tar -xvf libscl.tar
cd libscl/gpp
make
cd ..

# google test
wget https://github.com/google/googletest/archive/release-1.10.0.tar.gz
tar -xvf release-1.10.0.tar.gz
```

### Build

```sh
mkdir -p build
cd build
cmake ..
make
```

### Test

```sh
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -G "Unix Makefiles" 
make all
./test/jlsp_cpp_test 
```

## Usage

```sh
./jlst_cpp
```
