# Install

Requires UNIX environment

SRC

```shell
git clone git@github.com:MRCIEU/varGWAS.git
cd varGWAS
```

Load compiler (optional). Tested with GCC v5 & v6.

```shell
# BC4
module load languages/gcc/9.3.0
module load tools/cmake/3.20.0
module load OpenBLAS/0.2.20-GCC-6.4.0-2.28  
```

Libraries

- Assumes you already have [armadillo](http://arma.sourceforge.net/download.html) installed and dependencies

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
