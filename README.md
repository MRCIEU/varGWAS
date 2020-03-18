# JLST C++

Implementation of the JLST method (Staley et al) using C++

```
Staley, J. R., Windmeijer, F., Suderman, M., Smith, G. D., & Tilling, K. (n.d.). A robust mean and variance test with application to epigenome-wide association studies. https://doi.org/10.1101/2020.02.06.926584
```

# Install

## SRC

```sh
git clone git@ieugit-scmv-d0.epi.bris.ac.uk:ml18692/jlst_cpp.git
cd jlst_cpp
```

## Libraries

```sh
mkdir -p includes
cd includes
wget http://bitbucket.org/gavinband/bgen/get/master.tar.gz
tar -xvf master.tar.gz
cd ..
```

## Build

```sh
mkdir -p build
cd build
cmake ..
make
```

# Usage

```sh
./jlst_cpp ../includes/bgen/example/example.v11.bgen default
```
