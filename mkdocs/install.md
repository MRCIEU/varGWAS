# Install

## Precompiled binary

The precompiled binary for Linux systems can be downloaded from [GitHub](https://github.com/MRCIEU/varGWAS/releases). This is the simplest method and will work for most users.

## Build from source

Obtain source

```shell
git clone git@github.com:MRCIEU/varGWAS.git
cd varGWAS
```

Load compiler (may be necessary on HPC systems). Tested with GCC v7 & v9.

```shell
# BlueCrystal Phase 4
module load languages/gcc/9.3.0
module load tools/cmake/3.20.0
```

Build dependencies

```shell
bash lib.sh
```

Configure cmake

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
```

Build

```shell
make
```

Run

```shell
./bin/varGWAS
```

## Docker

Image

```shell
# pull image from Dockerhub
docker pull mrcieu/vargwas
### OR ###
# Build image from source
docker build -t vargwas .
```

Run

```shell
docker run \
-it \
-v $PWD:/home \
mrcieu/vargwas
```
