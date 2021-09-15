FROM gcc:6.5

# install cmake
RUN wget https://github.com/Kitware/CMake/releases/download/v3.19.6/cmake-3.19.6-Linux-x86_64.tar.gz && \
    tar -xvf cmake-3.19.6-Linux-x86_64.tar.gz
ENV PATH=$PATH:/cmake-3.19.6-Linux-x86_64/bin
ENV CMAKE_ROOT=/cmake-3.19.6-Linux-x86_64

# copy repo to container
COPY . /app
WORKDIR /app

# install libs
RUN bash lib.sh

# build executable
RUN mkdir /app/build && cd /app/build && cmake .. -DCMAKE_BUILD_TYPE=Release && make

# launch app
ENTRYPOINT ["/app/build/bin/varGWAS"]
