FROM ubuntu:20.04
ENV TZ=Europe/London
ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /tmp
RUN mkdir ~/apps

# install packages
RUN apt update && \
    apt install -y \
    build-essential \
    cmake \
    git \
    wget \
    r-base \
    python3

# alias python
RUN ln -s /usr/bin/python3 /usr/bin/python

# install R-packages
RUN R -e "install.packages(c('quadprog', 'pwr', 'cqrReg', 'data.table', 'optparse', 'quantreg', 'BEDMatrix', 'broom', 'genpwr', 'dplyr', 'remotes'))"
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
RUN ./waf configure --prefix=/usr
RUN ./waf
RUN ./waf install
RUN ln -s /usr/bin/qctool_v2.2.0 /usr/bin/qctool

# install BGENIX
WORKDIR /tmp
RUN wget http://code.enkre.net/bgen/tarball/release/bgen.tgz
RUN tar -xvf bgen.tgz
WORKDIR /tmp/bgen.tgz
RUN ./waf configure --prefix=/usr
RUN ./waf
RUN ./waf install

# install varGWAS
WORKDIR /usr/bin
RUN wget https://github.com/MRCIEU/varGWAS/releases/download/1.2.3/varGWAS
RUN chmod a+x ./varGWAS

# install QUAIL
WORKDIR /tmp
RUN git clone https://github.com/qlu-lab/QUAIL
RUN mv QUAIL ~/apps

# set working env
WORKDIR /app