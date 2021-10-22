# Install

Requires UNIX environment

SRC

```shell
git clone git@github.com:MRCIEU/varGWAS.git
cd varGWAS
```

Load compiler (optional). Tested with GCC v7 & v9.

```shell
# BC4
module load build/gcc-7.2.0
module load tools/cmake/3.20.0
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
CC=/mnt/storage/software/languages/gcc-7.2.0/bin/gcc \
CXX=/mnt/storage/software/languages/gcc-7.2.0/bin/g++ \
cmake .. -DCMAKE_BUILD_TYPE=Release



```

# Docker

Build image

```shell
docker build -t vargwas .
```
