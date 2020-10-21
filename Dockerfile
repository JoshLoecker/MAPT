# install guppy with the following version number
FROM alpine:3.6 as guppy_installer
MAINTAINER joshua.loecker@usda.gov

ARG VERSION_NUMBER=4.0.14

# specify and open a working directory
ARG INSTALL_DIR=/guppy
WORKDIR $INSTALL_DIR

# download and install guppy
RUN apk update && \
    apk add wget && \
    echo "" && \
    echo Downloading ont-guppy && \
    # switch the comment on the two lines below to use GPU instead of CPU
    wget --quiet --show-progress --progress=bar:force --no-check-certificate https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_${VERSION_NUMBER}_linux64.tar.gz  && \
    # wget --quiet --show-progress --no-check-certificate https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_${VERSION_NUMBER}_linux64.tar.gz && \
    tar -xf ont-guppy*.tar.gz && \
    rm ont-guppy*.tar.gz && \
    rm -rf /var/cache/apk/*

# get our conda environment running
FROM continuumio/miniconda3
MAINTAINER joshua.loecker@usda.gov

# set a working directory
WORKDIR /pipeline

# create a volume to store data
VOLUME ["/results/", "/data/", "/alignment_file.fasta"]

# copy the guppy downloader from previous workspace
# copy our environment file into the container
COPY --from=guppy_installer /guppy/ /pipeline
COPY environment.yml .

# we are doing a series of things here to reduce the final image size that is generated

# 2) create our pipeline environment, and update conda
# 3)
# 4)
#   - this includes the Snakefile, configuration file, environment scripts, and python scripts
# 5) remove the github folder as it is no longer needed
# 6) remove files that were loaded during the  `apt update && apt install` to reduce image size
# 7) create a link so we are able to redirect STDOUT from logs onto the screen

# update our package list and install a package needed for guppy, and git to clone the snakemake repo
RUN apt update && \
    apt install --yes git libidn11 && \
    git clone --branch master https://github.com/JoshLoecker/pipeline && \
    rm -rf /var/lib/apt/lists/* && \
    ln -s /pipeline/ont-guppy*/bin/guppy_aligner /usr/local/bin/guppy_aligner && \
    ln -s /pipeline/ont-guppy*/bin/guppy_barcoder /usr/local/bin/guppy_barcoder && \
    ln -s /pipeline/ont-guppy*/bin/guppy_basecall_server /usr/local/bin/guppy_basecall_server && \
    ln -s /pipeline/ont-guppy*/bin/guppy_basecaller /usr/local/bin/guppy_basecaller && \
    ln -s /pipeline/ont-guppy*/bin/guppy_basecaller_1d2 /usr/local/bin/guppy_basecaller_1d2 && \
    ln -s /pipeline/ont-guppy*/bin/guppy_basecaller_supervisor /usr/local/bin/guppy_basecaller_supervisor && \
    ln -sf  /dev/stdout /var/ && \
    mv /pipeline/pipeline/Snakefile /pipeline && \
    mv /pipeline/pipeline/config.yml /pipeline && \
    mv /pipeline/pipeline/envs /pipeline/envs/ && \
    mv /pipeline/pipeline/scripts /pipeline/scripts/ && \
    rm -r /pipeline/pipeline/ && \
    conda env create -f environment.yml && \
    conda update -n base -c defaults conda

# activate our new environment
SHELL ["conda", "run", "-n", "pipeline", "/bin/bash", "-c"]

# run the code contained in 'run_pipeline.py` when the `pipeline` environment is started
ENTRYPOINT ["conda", "run", "-n", "pipeline", "snakemake", "-j", "all"]
