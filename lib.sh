# script to build libraries
set -euxo pipefail
LIB_DIR=$(pwd)/lib
mkdir -p "$LIB_DIR"

# bgen
cd "$LIB_DIR"
curl -L http://code.enkre.net/bgen/tarball/release/bgen.tgz >bgen.tar.gz
tar -xf bgen.tar.gz
mv bgen.tgz bgen
cd bgen
./waf configure
./waf

# csv reader
# TODO replace with Boost
cd "$LIB_DIR"
curl -L https://github.com/ben-strasser/fast-cpp-csv-parser/archive/master.zip >fast-cpp-csv-parser.zip
unzip fast-cpp-csv-parser.zip
mv fast-cpp-csv-parser-master fast-cpp-csv-parser

# eigen (requires >=3.4)
cd "$LIB_DIR"
curl -L https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz >eigen-3.4.0.tar.gz
tar -xf eigen-3.4.0.tar.gz

# cxxopts
cd "$LIB_DIR"
curl -L https://github.com/jarro2783/cxxopts/archive/v2.2.0.tar.gz >cxxopts-2.2.0.tar.gz
tar -xf cxxopts-2.2.0.tar.gz

# spdlog
cd "$LIB_DIR"
git clone https://github.com/gabime/spdlog.git
cd spdlog && mkdir build && cd build
cmake .. && make -j

# google test
cd "$LIB_DIR"
curl -L https://github.com/google/googletest/archive/release-1.10.0.tar.gz >release-1.10.0.tar.gz
tar -xf release-1.10.0.tar.gz