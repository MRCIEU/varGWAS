FROM gcc:8.4

# install cmake
RUN apt-get update && apt-get -y install cmake protobuf-compiler

# copy repo to container
COPY . /app
WORKDIR /app

# install libs
RUN bash lib.sh

# build executable
RUN mkdir -p /app/build
WORKDIR /app/build
RUN cmake .. -DCMAKE_BUILD_TYPE=Release
RUN make

# launch app
WORKDIR /app/build/bin/
CMD ["jlst_cpp"]