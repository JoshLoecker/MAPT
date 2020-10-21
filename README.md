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
`data_files` folder. These folders can be named as you please, but this guide 
will use these respective names
It is important to note that under the `data_files` folder, a folder **must** 
be named `fast5`. If this is not done, guppy_basecaller will not be able to pick up **any** `.fast5` files, even if they are present, as it looks for `.fast5` files under the `fast5` folder. As such, your `data_files` folder structure may look as follows:
<br>
```
home
| -- Rob
    | -- Projects
        | -- run_1
            | -- data_files
                | -- fast5
                    | -- file_1.fast5
                    | -- file_2.fast5
                    | -- file_3.fast5
                | -- silva_alignment_file.fasta
                | -- some_other_file.txt
```
The `results` folder does not need to exist, but it can if you would like. If 
it does not exist, It will be created by Singularity/Docker when starting the 
container
### <ins>Singularty</ins>
1. If you are working on a cluster, Singularity may already be installed; this 
is the case with SciNet. If Singularity is not installed, follow the 
instructions: [Install Singularity](https://singularity.lbl.gov/install-linux)
2. If you would like to pull the CPU version of the pipeline, executing the 
following command in a terminal window:<br>
	`singularity pull docker://joshloecker/pipeline_cpu:latest`
3. The GPU version of the pipeline can be obtained by executing this command in
 a terminal window:<br>
	`singularity pull docker://joshloecker/pipeline_gpu:latest`
4. Then, set up a few environment variables:
    ```
    results="/path/to/results"
    data="/path/to/data"
    alignment_name="name_of_alignment_file.fasta"
    ```
5. To run the container, execute the following command:
    ```
    singularity run \
    --bind "${results}":/results/ \
    --bind "${data}":/data_files/ \
    --bind "${data}/${alignment_name}:/alignment_file.fasta \
    pipeline:latest
	```
### <ins>Docker</ins>
1. If you chose to work with Docker, ensure Docker is already installed on your 
system. If it is not, follow the instructions: 
[Install Docker](https://docs.docker.com/get-docker/)
2. Pull the container by executing the following command in a terminal 
window:<br>
	`docker pull joshloecker/pipeline:latest`
3. To start, set up a few variables in the terminal
    ```
    results="/path/to/results"
    data="/path/to/data"
    alignment_name="name_of_alignment_file.fasta"
   ```

4. Then, run the container using the following command. This can safely be copied
and pasted, assuming step 3 has been done
    ```
    docker run \
    --mount 'type=bind,source="${results}",target=/results/' \
    --mount 'type=bind,source="${data}":/data/,readonly' \
    --mount 'type=bind,source="${data}/${alignment_name}",target=/alignment_file.fasta,readonly' \
    joshloecker/pipeline_cpu:latest
    ```
   If you are using the GPU version of the pipeline, the final line 
   (`joshloecker/pipeline_cpu:latest`) should be changed to 
   `joshloecker/pipeline_gpu:latest`

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
it may be ran again.


## Known Errors
Deleting files from the `.temp` folder will cause the pipeline to regenerate 
these files, along with any output downstream.
If you experience any errors, contact joshua.loecker@usda.gov for assistance

<!--stackedit_data:
eyJoaXN0b3J5IjpbOTIyMDg4NzE1XX0=
-->
