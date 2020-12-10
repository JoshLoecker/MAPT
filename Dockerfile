FROM continuumio/miniconda3
MAINTAINER joshua.loecker@usda.gov
ARG VERSION_NUMBER=4.2.2

ENV alignment_name=${alignment_name}
ENV basecall_callers=${basecall_callers}
ENV threads_per_caller=${threads_per_caller}
VOLUME /results /data
# ADD https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_${VERSION_NUMBER}_linux64.tar.gz /guppy.tar.gz

# Remove after docker is good to go
COPY dockerBuild/ont-guppy-cpu guppy
COPY environment.yml environment.yml
COPY dockerBuild/pipeline /data/

# unpack guppy, link it to docker environment bash so it can be called anywhere
# add the these three lines when docker is running
# tar -xf guppy.tar.gz && \
    # rm guppy.tar.gz && \
    # mv /ont-guppy-cpu /guppy && \
RUN ln -s /guppy/bin/guppy_aligner /usr/local/bin/guppy_aligner && \
    ln -s /guppy/bin/guppy_barcoder /usr/local/bin/guppy_barcoder && \
    ln -s /guppy/bin/guppy_basecall_server /usr/local/bin/guppy_basecall_server && \
    ln -s /guppy/bin/guppy_basecaller /usr/local/bin/guppy_basecaller && \
    ln -s /guppy/bin/guppy_basecaller_1d2 /usr/local/bin/guppy_basecaller_1d2 && \
    ln -s /guppy/bin/guppy_basecaller_supervisor /usr/local/bin/guppy_basecaller_supervisor


# clone github repo when docker is running
#git clone --branch master https://github.com/JoshLoecker/pipeline && \
RUN apt update && \
    apt --yes --no-install-recommends install g++ git libidn11 m4 autoconf automake libtool && \
    # update conda environment
    #conda update -n base -c defaults conda && \
    # create conda environment with our environment file
    #conda env create --file /pipeline/environment.yml && \
    conda env create --file environment.yml && \
    # remove unneeded apt files to reduce image size
    apt --yes purge g++ git && \
    apt --yes autoremove && \
    rm -rf /var/lib/apt/lists/*


ENV PATH /opt/conda/envs/pipeline/bin:$PATH
RUN /bin/bash -c  "source activate pipeline"
CMD ["pip", "install", "networkx", "IsoCon", "isONclust"]
# RUN conda activate pipeline && \
#     conda --version
#     #pip install networkx==2.3 IsoCon==0.3.* isONclust==0.*

# # start our conda environment `pipeline`, and call `snakemake`
ENTRYPOINT cd /data && conda run -n pipeline snakemake -j all
