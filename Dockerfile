FROM ubuntu:16.04

RUN apt update && \
    apt install -y git wget make g++ zlib1g-dev python-pip python3-pip && \
    rm -rf /var/lib/apt/lists/*

# Install bedtools v2.27.1
RUN git clone --branch v2.27.1 --single-branch https://github.com/arq5x/bedtools2.git && cd bedtools2 && make && make install && cd ../ && rm -rf bedtools2

RUN wget -q https://github.com/bedops/bedops/releases/download/v2.4.32/bedops_linux_x86_64-v2.4.32.tar.bz2 && mkdir -p bedops_tmp && tar -C bedops_tmp -jxf bedops_linux_x86_64-v2.4.32.tar.bz2 && mv bedops_tmp/bin/* /usr/local/bin && rm -rf bedops_tmp bedops_linux_x86_64-v2.4.32.tar.bz2

RUN pip3 install --user six requests python-dateutil

RUN mkdir -p /app/
#COPY ./src /app/src
#COPY ./input /app/input

#ENTRYPOINT ["/bin/bash","-c"]