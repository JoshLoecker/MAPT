Pipeline (MAPT)
===============

Naming
------
I would like to name this project as it has been my first, and most in-depth, project to date. 
Due to the nature of this project revolving around mapping noisy microbial reads (hopefully) down to the species level, 
the following name seemed appropriate: Microbial Automated Processing Tool, or MAPT for short.


Usage
-----

This is a workflow/pipeline aimed to streamline post-processing of noisy reads from [Oxford Nanopore MinION Sequencing Data](https://www.nanoporetech.com). The following software will need to be installed in a conda environment. A startup-script is provided in the GitHub repository. This is explained further on this page.

The following software will be used (not necessarily in this order):
- Guppy Basecalling
- Guppy Barcoding
- Cutadapt
- NanoFilt
- isONclust
- spoa
- NanoPlot
- Plotly Graphs

The page to the GitHub and Docker repositories are as follows, if they are needed.
- GitHub: [https://github.com/JoshLoecker/pipeline](https://github.com/JoshLoecker/pipeline)
- Docker: [https://hub.docker.com/repository/docker/joshloecker/pipeline](https://hub.docker.com/repository/docker/joshloecker/pipeline)

Installation & Setup
------------

This project was built and testing on the following versions of Singularity and Docker. Reproducibility and stability cannot be guaranteed on earlier versions.
1. Singularity (version 3.6.4-1.el7 or higher)
2. Docker (version 20.10.0 or higher)

[Download singularity here](https://singularity.lbl.gov/)  
[Download docker here](https://www.docker.com/products/docker-desktop)

Git and Conda are required to download the pipeline and create a new Conda environment with the required software

1. Start by cloning the GitHub repo  
    a. `git clone https://github.com/JoshLoecker/pipeline`
    b. This will create a new folder `pipeline` in your current directory. Move this to the folder you would like it to be contained in
2. Within the `pipeline` folder, edit the `environment.yaml` file.  
    a. Change the last line to the location where you would like to store your conda environments  
    b. On SciNet, storage space in your home directory is limited. It is [strongly recommended](https://scinet.usda.gov/guide/ceres/#quotas-on-home-and-project-directories) to store environments in your `/project` directory.  
    c. By saving environments in your `/project` directory, it is possible for everyone on your team to use the same environment, even if multiple people are running jobs at the same time.
    d. You may also change the name of the pipeline by changing `name: mapt_pipeline` to `name: MY_NAME` if you choose
   
If this is the first time setting up the conda environment for your group, continue to Step 3
If step 2 has been completed for your group already, you are ready to activate the conda environment. Proceed to Step 4

3. Navigate to the `setup` folder where you first downloaded the pipeline (from Step 1)  
    a. Call the `setup.sh` script by running `./setup.sh`  
    b. This will guide you through a simple setup to install the conda environment  
4. If you know where the `prefix: ` to the conda environment was set, simply type `source activate /path/to/conda/env/name`  
    a. To deactivate the environment, type `conda deactivate`.  
5. To be able to use the Guppy suite of tools, the 
5. The final step is to edit several lines within the `pipeline/config.yaml` file  
    a. First, set the `results`, `data`, `reference_database`, and `guppy_container` to their appropriate locations  
	1) `results` is where you would like results of the pipeline to be stored  
	2) `data` holds your fast5/fastq files for Guppy  
	3) `reference_database` is the database you will be using with MiniMap for alignments  
	4) `guppy_container` is the location of the guppy container you will be using  
	
	b. Then set any other values required under the `DEFAULT VALUES` section. If these are not changed, they will remain as-is during the pipeline run



Notes to Future Maintainers
---------------------------
1. Singularity and Docker must be installed on the same machine to update guppy
2. Building the Guppy singularity image was first done by building a docker container  
	a. `docker build --tag [YOUR TAG] .`  
    b. This was done simply because I was most familiar with docker containers  
    c. It may be smart to move the Dockerfile in the `pipeline` repository to a Singularity file  
3. The singularity container is built in the following manner  
	a. `singularity build --sandbox docker-daemon://[YOUR TAG FROM STEP 2]`  
	b. The `docker-daemon` is used for a local docker image. Local images are generally preferred. This means we do not have to upload the resulting container to Dockerhub, then download it to our local machine  
	1) `singularity build --sandbox docker://[YOUR TAG]` will download a docker container from Dockerhub, if this is preferred.  
    
	c. 
