FROM continuumio/miniconda3
MAINTAINER joshua.loecker@usda.gov
ARG VERSION_NUMBER=4.2.2
ENV alignment_name=${alignment_name}
VOLUME /results /data
# ADD https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_${VERSION_NUMBER}_linux64.tar.gz /guppy.tar.gz
COPY guppy.tar.gz guppy.tar.gz

# unpack guppy, link it to docker environment bash so it can be called anywhere
RUN tar -xf guppy.tar.gz && \
    rm guppy.tar.gz && \
    mv /ont-guppy-cpu /guppy && \
    ln -s /guppy/bin/guppy_aligner /usr/local/bin/guppy_aligner && \
    ln -s /guppy/bin/guppy_barcoder /usr/local/bin/guppy_barcoder && \
    ln -s /guppy/bin/guppy_basecall_server /usr/local/bin/guppy_basecall_server && \
    ln -s /guppy/bin/guppy_basecaller /usr/local/bin/guppy_basecaller && \
    ln -s /guppy/bin/guppy_basecaller_1d2 /usr/local/bin/guppy_basecaller_1d2 && \
    ln -s /guppy/bin/guppy_basecaller_supervisor /usr/local/bin/guppy_basecaller_supervisor

# download GitHub repo
RUN apt update && \
    apt --yes --no-install-recommends install g++ git libidn11 m4 autoconf automake libtool && \
    # get github repo to the folder
    git clone --branch master https://github.com/JoshLoecker/pipeline && \
    # update conda environment
    conda update -n base -c defaults conda && \
    # create conda environment with our environment file
    ls / && \
    conda env create --file /pipeline/environment.yml && \
    # remove installed packages, as they are no longer needed
    apt --yes purge g++ git && \
    # remove unneeded apt files to reduce image size
    rm -rf /var/lib/apt/lists/*

RUN conda activate pipeline && \
    pip install

# # start our conda environment `pipeline`, and call `snakemake`
ENTRYPOINT cd pipeline && conda run -n pipeline snakemake -j all
