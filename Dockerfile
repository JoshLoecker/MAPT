FROM continuumio/miniconda3
MAINTAINER joshua.loecker@usda.gov
ARG VERSION_NUMBER=4.2.2

ENV alignment_path=""
ENV num_basecallers=""
ENV threads_per_caller=""
ENV basecall_configuration=""
ENV barcode_kit=""
ENV cutadapt_trim_error_rate=0.15
ENV cutadapt_trim_three_prime_adapter=""
ENV cutadapt_trim_five_prime_adapter=""
ENV nanofilt_filtering_min=1000
ENV nanofilt_filtering_max=2000
ENV basecall=""

# Copy files
COPY dockerBuild/ont-guppy /guppy
COPY Dockerfile /workflow/
COPY Snakefile /workflow/
COPY config.yaml /workflow/
COPY environment.yml /workflow/
COPY scripts /workflow/scripts

# unpack guppy, link it to docker environment bash so it can be called anywhere
# add the these three lines when docker is running
# tar -xf guppy.tar.gz && \
    # rm guppy.tar.gz && \
    # mv /ont-guppy-cpu /guppy && \

# clone github repo when docker is running
#git clone --branch master https://github.com/JoshLoecker/pipeline workflow && \
RUN conda env create --file /workflow/environment.yml && \
    #conda update -n base -c defaults conda && \
    # create conda environment with our environment file
    # link guppy to docker environment so it can be called anywhere
    ln -s /guppy/bin/guppy_aligner /usr/local/bin/guppy_aligner && \
    ln -s /guppy/bin/guppy_barcoder /usr/local/bin/guppy_barcoder && \
    ln -s /guppy/bin/guppy_basecall_server /usr/local/bin/guppy_basecall_server && \
    ln -s /guppy/bin/guppy_basecaller /usr/local/bin/guppy_basecaller && \
    ln -s /guppy/bin/guppy_basecaller_supervisor /usr/local/bin/guppy_basecaller_supervisor && \
    ln -s /workflow/data /data && \
    ln -s /workflow/results /results
ENV PATH /opt/conda/envs/pipeline/bin:$PATH

RUN apt update && \
    # install required packages
    apt --yes --no-install-recommends install g++ libidn11 m4 autoconf automake dh-autoreconf && \
    # clone parasail (required by IsoCon)
    git clone https://github.com/jeffdaily/parasail-python /parasail && \
    # remove unneeded apt files to reduce image size
    apt --yes purge git && \
    apt --yes autoremove && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /tmp/* && \
    rm -rf /var/tmp/*

# activate our conda environment
RUN /bin/bash -c  "source activate pipeline" && \
    # make the alignment_name, num_basecallers, and threads_per_caller available for python as environment variables
    # build parasail (required by IsoCon, could not get this to work with pip install)
    cd /parasail && \
    python setup.py bdist_wheel && \
    pip install IsoCon isONclust

# start our conda environment `pipeline`, and call `snakemake`
# adding `-n` or `--dry-run` after calling docker run joshloecker/pipeline will perform a dry-run of the pipeline
ENTRYPOINT ["sh", "-c", "conda", "run", "--name", "pipeline", "--cwd", "/workflow", "snakemake", "-j", "all", "--config", "basecall", "=", "$basecall"]
VOLUME /results
VOLUME /data
