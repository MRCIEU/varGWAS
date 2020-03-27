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

# csv reader 
curl -L https://github.com/ben-strasser/fast-cpp-csv-parser/archive/713c5fd.zip > fast-cpp-csv-parser-713c5fd2ba1b6d145296a21fc7f9dee576daaa4f.zip
unzip fast-cpp-csv-parser-713c5fd2ba1b6d145296a21fc7f9dee576daaa4f.zip

# eigen
curl https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz > eigen-3.3.7.tar.gz
tar -xvf eigen-3.3.7.tar.gz 

# cxxopts
curl -L https://github.com/jarro2783/cxxopts/archive/v2.2.0.tar.gz > cxxopts-2.2.0.tar.gz
tar -xvf cxxopts-2.2.0.tar.gz

# google logging
curl -L https://github.com/google/glog/archive/v0.4.0.tar.gz > glog-0.4.0.tar.gz
tar -xvf glog-0.4.0.tar.gz
cd glog-0.4.0
mkdir bin
cmake ..
make

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
