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

Installation
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
5. The final step is to 
