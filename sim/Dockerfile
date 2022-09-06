FROM ubuntu:22.04
ENV TZ=Europe/London
ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /tmp
RUN mkdir ~/apps

# install packages
RUN apt update && \
    apt install -y \
    build-essential \
    r-base \
    python3

# install R-packages
RUN R -e "install.packages(c('pwr', 'cqrReg', 'data.table', 'optparse', 'quantreg', 'BEDMatrix', 'broom', 'genpwr', 'dplyr', 'remotes'))"
RUN R -e "remotes::install_github('jrs95/jlst')"
RUN R -e "remotes::install_github('mrcieu/varGWASR')"

# install OSCA
RUN wget https://yanglab.westlake.edu.cn/software/osca/download/osca-0.46.1-linux-x86_64.zip
RUN unzip osca-0.46.1-linux-x86_64.zip
RUN mv osca-0.46.1-linux-x86_64 /usr/bin/osca

# install DRM
RUN git clone http://github.com/drewmard/DRM
RUN mv DRM ~/apps

# install qctool
RUN wget https://code.enkre.net/qctool/zip/release/qctool.tgz
RUN unzip qctool.tgz
WORKDIR /tmp/qctool
RUN ./waf configure
RUN ./waf

# install BGENIX
WORKDIR /tmp
RUN wget http://code.enkre.net/bgen/tarball/release/bgen.tgz
RUN tar -xvf bgen.tgz
WORKDIR /tmp/bgen
RUN ./waf configure
RUN ./waf

# install varGWAS
WORKDIR /tmp
RUN git clone git@github.com:MRCIEU/varGWAS.git
WORKDIR /tmp/varGWAS
RUN bash lib.sh
RUN mkdir -p build
WORKDIR /tmp/varGWAS/build
RUN cmake .. -DCMAKE_BUILD_TYPE=Release

# install QUAIL
WORKDIR /tmp
RUN git clone https://github.com/qlu-lab/QUAIL


WORKDIR /app