# script to build libraries
set -euxo pipefail
LIB_DIR=$(pwd)/lib
mkdir -p "$LIB_DIR"

# bgen
cd "$LIB_DIR"
curl -L http://code.enkre.net/bgen/tarball/release/bgen.tgz >bgen.tar.gz
tar -xvf bgen.tar.gz
mv bgen.tgz bgen
cd bgen
./waf configure
./waf

# boost
cd "$LIB_DIR"
curl -L https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.bz2 >boost_1_75_0.tar.bz2
tar --bzip2 -xf boost_1_75_0.tar.bz2

# zstd
cd "$LIB_DIR"
curl -L https://github.com/facebook/zstd/archive/v1.1.0.tar.gz > zstd-1.1.0.tar.gz
tar -xvf zstd-1.1.0.tar.gz
cd zstd-1.1.0/build/cmake
mkdir build
cd build
cmake ..
make

# quantile regression
cd "$LIB_DIR"
curl -L http://www.aronaldg.org/webfiles/libscl/libscl.tar >libscl.tar
tar -xvf libscl.tar
cd libscl/gpp
make

# csv reader
cd "$LIB_DIR"
curl -L https://github.com/ben-strasser/fast-cpp-csv-parser/archive/master.zip >fast-cpp-csv-parser.zip
unzip fast-cpp-csv-parser.zip
mv fast-cpp-csv-parser-master fast-cpp-csv-parser

# eigen (requires >=3.4)
cd "$LIB_DIR"
curl -L https://gitlab.com/libeigen/eigen/-/archive/master/eigen-master.tar.gz >eigen-master.tar.gz
tar -xvf eigen-master.tar.gz

# cxxopts
cd "$LIB_DIR"
curl -L https://github.com/jarro2783/cxxopts/archive/v2.2.0.tar.gz >cxxopts-2.2.0.tar.gz
tar -xvf cxxopts-2.2.0.tar.gz

# spdlog
cd "$LIB_DIR"
git clone https://github.com/gabime/spdlog.git
cd spdlog && mkdir build && cd build
cmake .. && make -j

# google test
cd "$LIB_DIR"
curl -L https://github.com/google/googletest/archive/release-1.10.0.tar.gz >release-1.10.0.tar.gz
tar -xvf release-1.10.0.tar.gz

# ThreadPool
cd "$LIB_DIR"
curl -L https://github.com/progschj/ThreadPool/archive/master.zip >ThreadPool.zip
unzip ThreadPool.zip
mv ThreadPool-master ThreadPool