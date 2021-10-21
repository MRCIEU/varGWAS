# Install

Requires UNIX environment

SRC

```shell
git clone git@github.com:MRCIEU/varGWAS.git
cd varGWAS
```

Load compiler (optional). Tested with GCC v9.

```shell
# BC4
module load languages/gcc/9.3.0
module load tools/cmake/3.20.0
module load ScaLAPACK/2.0.2-gompic-2016.10-OpenBLAS-0.2.19-LAPACK-3.6.1
module load HDF5/1.8.17-foss-2016b
```

Libraries

```shell
bash lib.sh
```

Build

```shell
mkdir -p build
cd build

# use default compiler path
cmake .. -DCMAKE_BUILD_TYPE=Release

### OR ###

# use custom compiler path
CC=/mnt/storage/software/languages/gcc-9.3/bin/gcc \
CXX=/mnt/storage/software/languages/gcc-9.3/bin/g++ \
cmake .. -DCMAKE_BUILD_TYPE=Release

make
```

# Docker

Build image

```shell
docker build -t vargwas .
```
