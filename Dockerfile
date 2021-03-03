FROM gcc:5.5

# install libs
#RUN apt-get update && apt-get -y install zlib1g-dev

# install cmake binaries
RUN wget https://github.com/Kitware/CMake/releases/download/v3.19.6/cmake-3.19.6-Linux-x86_64.tar.gz && \
    tar -xvf cmake-3.19.6-Linux-x86_64.tar.gz
ENV PATH=/cmake-3.19.6-Linux-x86_64/bin:$PATH
ENV CMAKE_ROOT=/cmake-3.19.6-Linux-x86_64

# copy repo to container
COPY . /app
WORKDIR /app

# install libs
RUN bash lib.sh

# build executable
RUN mkdir -p /app/build
WORKDIR /app/build
RUN cmake .. -DCMAKE_BUILD_TYPE=Release && make

# launch app
WORKDIR /app/build/bin
CMD ["jlst_cpp"]