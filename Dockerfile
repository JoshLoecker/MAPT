# install guppy with the following version number
FROM alpine:3.6 as guppy_downloader
MAINTAINER joshua.loecker@usda.gov

ARG VERSION_NUMBER=4.2.2

# specify and open a working directory
ARG INSTALL_DIR=/guppy

#ADD https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_${VERSION_NUMBER}_linux64.tar.gz /guppy
# download and install guppy
RUN mkdir ${INSTALL_DIR} && \
    cd ${INSTALL_DIR} && \
    apk update && \
    apk add wget && \
    echo "" && \
    echo Downloading ont-guppy && \
    # switch the comment on the two lines below to use GPU instead of CPU
    # curl -L -o guppy.tar.gz https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_${VERSION_NUMBER}_linux64.tar.gz && \
    wget --quiet --show-progress --progress=bar:force --no-check-certificate -O guppy.tar.gz https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_${VERSION_NUMBER}_linux64.tar.gz && \
    # wget --quiet --show-progress --no-check-certificate https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_${VERSION_NUMBER}_linux64.tar.gz && \
    echo Done && \
    tar -xf guppy.tar.gz && \
    rm guppy.tar.gz && \
    apk del wget && \
    rm -rf /var/cache/apk/*

# get our conda environment running
FROM continuumio/miniconda3
MAINTAINER joshua.loecker@usda.gov

# set a working directory
WORKDIR /pipeline
ENV alignment_name=${alignment_name}

# create a volume to store data
VOLUME ["/results/", "/data/", "/alignment_file.fasta"]

# copy the guppy downloader from previous workspace
# copy our environment file into the container
COPY --from=guppy_downloader /guppy/ /pipeline
COPY environment.yml .

# update our package list and install a package needed for guppy, and git to clone the snakemake repo
RUN apt update && \
    apt install --yes --no-install-recommends g++ git libidn11 && \
    git clone --branch master https://github.com/JoshLoecker/pipeline && \
    conda update -n base -c defaults conda && \
    conda env create -f environment.yml && \
    rm -rf /var/lib/apt/lists/* && \
    apt purge --yes g++ git

# move files from previous layer into this one for comand line use
RUN mv /pipeline/pipeline/Snakefile /pipeline && \
    mv /pipeline/pipeline/config.yml /pipeline && \
    mv /pipeline/pipeline/envs /pipeline/envs/ && \
    mv /pipeline/pipeline/scripts /pipeline/scripts/ && \
    rm -r /pipeline/pipeline/ && \
    # link guppy packages
    ln -s /pipeline/ont-guppy*/bin/guppy_aligner /usr/local/bin/guppy_aligner && \
    ln -s /pipeline/ont-guppy*/bin/guppy_barcoder /usr/local/bin/guppy_barcoder && \
    ln -s /pipeline/ont-guppy*/bin/guppy_basecall_server /usr/local/bin/guppy_basecall_server && \
    ln -s /pipeline/ont-guppy*/bin/guppy_basecaller /usr/local/bin/guppy_basecaller && \
    ln -s /pipeline/ont-guppy*/bin/guppy_basecaller_1d2 /usr/local/bin/guppy_basecaller_1d2 && \
    ln -s /pipeline/ont-guppy*/bin/guppy_basecaller_supervisor /usr/local/bin/guppy_basecaller_supervisor && \
    ln -sf  /dev/stdout /var/

# activate our new environment
SHELL ["conda", "run", "-n", "pipeline", "/bin/bash", "-c"]

# start our conda environment `pipeline`, and call `snakemake`
ENTRYPOINT ["conda", "run", "-n", "pipeline", "snakemake", "-j", "all"]
