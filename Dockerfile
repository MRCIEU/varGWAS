FROM ubuntu
ENV DEBIAN_FRONTEND=noninteractive

# install BLAS, lapack, cmake and compiler
RUN apt-get update && apt-get install -y build-essential cmake curl python git zlib1g-dev unzip libopenblas-dev liblapack-dev libarpack2-dev libsuperlu-dev

# copy repo to container
COPY . /app
WORKDIR /app

# install libs
RUN bash lib.sh

# build executable
RUN mkdir /app/build && cd /app/build && cmake .. -DCMAKE_BUILD_TYPE=Release && make

# launch app
ENTRYPOINT ["/app/build/bin/varGWAS"]