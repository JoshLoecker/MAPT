FROM ubuntu:latest as guppy
ARG version=4.4.1

# Ubuntu-20 GPU
ADD https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_${version}-1~focal_amd64.deb /guppy.deb

RUN apt update && \
    # install guppy
    apt install --yes /guppy.deb && \
    rm /guppy.deb && \
    # remove unnecessary items in image
    apt-get clean autoclean && \
    apt-get autoremove --yes && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
