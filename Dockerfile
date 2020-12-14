FROM continuumio/miniconda3
MAINTAINER joshua.loecker@usda.gov
ARG VERSION_NUMBER=4.2.2

ENV alignment_name=''
ENV num_basecallers=''
ENV threads_per_caller=''
# ADD https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_${VERSION_NUMBER}_linux64.tar.gz /guppy.tar.gz

# Remove after docker is good to go
COPY dockerBuild/ont-guppy-cpu /guppy
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
RUN ln -s /guppy/bin/guppy_aligner /usr/local/bin/guppy_aligner && \
    ln -s /guppy/bin/guppy_barcoder /usr/local/bin/guppy_barcoder && \
    ln -s /guppy/bin/guppy_basecall_server /usr/local/bin/guppy_basecall_server && \
    ln -s /guppy/bin/guppy_basecaller /usr/local/bin/guppy_basecaller && \
    ln -s /guppy/bin/guppy_basecaller_1d2 /usr/local/bin/guppy_basecaller_1d2 && \
    ln -s /guppy/bin/guppy_basecaller_supervisor /usr/local/bin/guppy_basecaller_supervisor


# clone github repo when docker is running
#git clone --branch master https://github.com/JoshLoecker/pipeline workflow && \
RUN apt update && \
    apt --yes --no-install-recommends install g++ git libidn11 m4 autoconf automake dh-autoreconf nano && \
    # update conda environment
    #conda update -n base -c defaults conda && \
    # create conda environment with our environment file
    #conda env create --file /workflow/environment.yml && \
    conda env create --file /workflow/environment.yml && \
    # clone parasail (required by IsoCon)
    git clone https://github.com/jeffdaily/parasail-python /parasail && \
    # remove unneeded apt files to reduce image size
    apt --yes purge git && \
    apt --yes autoremove && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /tmp/* && \
    rm -rf /var/tmp/*

ENV PATH /opt/conda/envs/pipeline/bin:$PATH

# activate our conda environment
RUN /bin/bash -c  "source activate pipeline" && \
    # make the alignment_name, num_basecallers, and threads_per_caller available for python as environment variables
    # build parasail (required by IsoCon, could not get this to work with pip install)
    cd /parasail && \
    python setup.py bdist_wheel && \
    pip install IsoCon isONclust

# start our conda environment `pipeline`, and call `snakemake`
# adding `-n` or `--dry-run` after calling docker run joshloecker/pipeline will perform a dry-run of the pipeline
CMD ["conda", "run", "--name", "pipeline", "--cwd", "/workflow", "snakemake", "-j", "all"]
VOLUME /results
VOLUME /data
