FROM continuumio/miniconda3
MAINTAINER joshua.loecker@usda.gov
ARG VERSION_NUMBER=4.2.2

# ADD https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_${VERSION_NUMBER}_linux64.tar.gz /guppy.tar.gz
COPY /dockerBuild/ont-guppy.tar.gz /guppy.tar.gz

RUN apt update && \
    # install required packages
    apt --yes --no-install-recommends install g++ libidn11 m4 autoconf automake dh-autoreconf && \
    # clone parasail (required by IsoCon) and the pipeline
    git clone https://github.com/jeffdaily/parasail-python /parasail && \
    git clone --branch master https://github.com/JoshLoecker/pipeline /workflow && \
    # build parasail (could not get this to work with pip install, had to go this way)
    cd /parasail && \
    python setup.py bdist_wheel && \
    cd / && \
    # set up conda environment
    conda env create --file /workflow/environment.yaml && \
    # unpack guppy, use symbolic link so it can be used anywhere
    tar -xf /guppy.tar.gz && \
    rm /guppy.tar.gz && \
    mv /ont-guppy* /guppy && \
    ln -s /guppy/bin/guppy_aligner /usr/local/bin/guppy_aligner && \
    ln -s /guppy/bin/guppy_barcoder /usr/local/bin/guppy_barcoder && \
    ln -s /guppy/bin/guppy_basecall_server /usr/local/bin/guppy_basecall_server && \
    ln -s /guppy/bin/guppy_basecaller /usr/local/bin/guppy_basecaller && \
    ln -s /guppy/bin/guppy_basecaller_supervisor /usr/local/bin/guppy_basecaller_supervisor && \
    # remove unneeded apt files to reduce image size
    apt --yes purge git && \
    apt --yes autoremove && \
    apt clean && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /tmp/* && \
    rm -rf /var/tmp/* && \
    # link conda so it can be used by singularity
    ln -s /opt/miniconda3/bin/conda /usr/bin/conda

# start our conda environment `pipeline`, and call `snakemake`
ENTRYPOINT ["conda", "run", "--name", "pipeline", "--cwd", "/workflow", "snakemake", "--cores", "all"]
VOLUME /results
VOLUME /data
