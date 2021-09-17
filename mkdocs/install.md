# Install

Requires UNIX environment

SRC

```shell
git clone git@github.com:MRCIEU/varGWAS.git
cd varGWAS
```

Libraries

```shell
bash lib.sh
```

Build

Tested with GCC v5 & v6. [Libscl](http://www.aronaldg.org/webfiles/libscl/) does not build with GCC >=v7.

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
module load tools/cmake/3.20.0
CC=/mnt/storage/software/languages/gcc-5.5.0/bin/gcc \
CXX=/mnt/storage/software/languages/gcc-5.5.0/bin/g++ \
cmake .. -DCMAKE_BUILD_TYPE=Release
```

# Docker

Build image

```shell
docker build -t varGWAS .
```