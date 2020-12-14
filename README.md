# Pipeline

## Naming
I would like to name this project as it has been my first, and most in-depth, 
project to date. Due to the nature of this project revolving around mapping 
noisy reads (hopefully) down to the species level, the following name seemed 
appropriate: Microbial Automated Processing Tool, or MAPT for short.

## Installation and Running
This project is meant to run using Singularity or Docker
### <ins>Folder Setup</ins>
Two paths are needed with this container: 1) A `results` folder, and 2) A 
`data` folder. These folders can be named as you please, but this guide 
will use these respective names.  

Under the `data` folder, a folder **must** 
be named `fast5`. If this is not done, guppy_basecaller will not be able to find
**any** `.fast5` files, even if they are present, as it looks for `.fast5`
files under the `fast5` folder.  

The `results` folder must exist before running the container. Your alignment file
should also be placed within the `data` folder.

Your `data` folder structure may look as follows:
```
home
| -- Rob
    | -- Projects
        | -- alignment_files
            | -- silva_alignment_file.fasta
            | -- another_alignment_file.fasta
        | -- run_1
            | -- data
                | -- fast5
                    | -- file_1.fast5
                    | -- file_2.fast5
                    | -- file_3.fast5
                | -- some_other_file.txt
            | -- results
```

### <ins>To Start</ins>
Set up a few variables in the terminal, using the above structure as an example 
Change this to values that make sense for your workflow

    results="/home/Rob/Projects/run_1/results"
    data="/home/Rob/Projects/run_1/data/"
    alignment_path="/home/Rob/Projects/alignment_files/silva_alignment_file.fasta"
    num_basecallers=3
    threads_per_caller=5
    perform_basecall=True

   The multiplication of `basecall_callers` and `threads_per_caller` should be
   very close to the number cores/threads you have on your machine, or reserved
   on a cluster.  
   If you do not have fast5 files (i.e. no need for basecalling),
   set the `perform_basecall` option to `False`.
  
   
### <ins>Singularty</ins>
1. If you chose to work with Singularity, ensure it is already installed on your 
system by running `singularity --version`. If `singularity version . . .` does not appear,
[install singularity here](https://singularity.lbl.gov/install-linux).  
Singularity is already installed on most clusters, such as SciNET

2. Next, pull the docker container  
    `singularity pull docker://joshloecker/pipeline:latest`

4. Assuming the `To Start` section was completed, the following command
can safely be copied and pasted. This will download the container, but it is
not yet running
    ```
    singularity create \
    --bind "${results}":/results/ \
    --bind "${data}":/data_files/ \
    --env alignment_path="${alignment_path}" \
    --env num_basecallers=${num_basecallers} \
    --env threads_per_caller=${threads_per_caller} \ 
    joshloecker/pipeline:latest
	```

5. To run the container, perform the following
   ```
   singularity run pipeline
   ```
   
6. If you would like to see a dry-run of the pipeline, append `--dry-run`
to the end of the `run` command. [Other Snakemake flags](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
can also be added in this manner.  
   ```
   singularity run joshloecker/pipeline --dry-run
   ```
    Note: Due to this workflow running in a container, not all flags have been tested
    nor are confirmed to work as expected.  
   
   
### <ins>Docker</ins>
1. If you chose to work with Docker, ensure Docker is already installed on your 
system by running `docker --version`. If `Docker version . . .` does not appear,
[install docker here](https://docs.docker.com/get-docker/)  

2. Then, run the container using the following command. This can safely be copied
and pasted, assuming the previous step has been completed
    ```
    docker create \
    --name=pipeline \
    --mount type=bind,source="${results}",target=/results \
    --mount type=bind,source="${data}",target=/data \
    -e alignment_path="${alignment_path}" \
    -e num_basecallers=${num_basecallers} \
    -e threads_per_caller=${threads_per_caller} \
    joshloecker/pipeline:latest
    ```
This will only download the container

3. To run the container, simply run
   ```
   docker run joshloecker/pipeline 
   ```
  
4. If you would like to see a dry-run of the pipeline, append `--dry-run`
to the end of the `run` command. [Other Snakemake flags](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
can also be added in this manner.  
   ```
   docker run joshloecker/pipeline --dry-run
   ```
    Note: Due to this workflow running in a container, not all flags have been tested
    nor are confirmed to work as expected. 

## Workflow
The following workflow will be completed, relatively in this order
1. [Oxford Nanopore Basecalling](https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revt_14dec2018) (i.e. guppy_basecaller)
2. [Oxford Nanopore Barcoding](https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revt_14dec2018) (i.e. guppy_barcoder)
3. [NanoQC](https://github.com/wdecoster/nanoQC)  
    a. NanoQC will be called independently on basecalling and barcoding results
4. [NanoPlot](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwiz89ql2cHsAhUjAp0JHUVoCFAQFjAAegQIARAC&url=https%3A%2F%2Fgithub.com%2Fwdecoster%2FNanoPlot&usg=AOvVaw00LEGNovoQzjS5KCUxwD0v)
5. [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
6. [NanoFilt](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwigp9S92cHsAhWSLc0KHYp8C-oQFjAAegQIBBAC&url=https%3A%2F%2Fgithub.com%2Fwdecoster%2Fnanofilt&usg=AOvVaw20npdGb-VRvmFH1SY6-P6C)
7. Oxford Nanopore Alignment (i.e. guppy_aligner)
8. [Minimap Aligner](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwjUtqXD2cHsAhVCGM0KHYsKAqcQFjAAegQIARAC&url=https%3A%2F%2Fgithub.com%2Flh3%2Fminimap2&usg=AOvVaw3UvK2vgr0fj_4GS68K8V26)
9. [VSearch Aligner](https://github.com/torognes/vsearch)
10. [Plotly](https://pypi.org/project/plotly/) Visuals

## Results
Results from the pipeline will be saved to the path you bound during the 
**Installation and Running** step (i.e. `/path/to/results/`). From here, you 
will be able to see each of the folders that are output, such as `basecalling`, 
`barcoding`, `visuals`, etc.

A `.temp` folder is generated under the `results` folder. This file contains
data that is not necessary to the general user. If this folder is deleted, the
next `snakemake` run will detect results in this folder is not present, and
the pipeline will be ran again.

## Known Errors
Deleting files from the `.temp` folder will cause the pipeline to regenerate 
these files, along with any output downstream.
If you experience any errors, contact [joshua.loecker@jacks.sdstate.edu](mailto:joshua.loecker@jacks.sdstate.edu) for assistance

<!--stackedit_data:
eyJoaXN0b3J5IjpbOTIyMDg4NzE1XX0=
-->
