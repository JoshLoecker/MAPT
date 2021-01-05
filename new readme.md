# Pipeline (MAPT)

### Naming
I would like to name this project as it has been my first, and most in-depth, project to date. 
Due to the nature of this project revolving around mapping noisy microbial reads (hopefully) down to the species level, 
the following name seemed appropriate: Microbial Automated Processing Tool, or MAPT for short.


### Usage

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

### Installation

This project was built and testing on the following versions of Singularity and Docker. Reproducibility and stability cannot be guaranteed on earlier versions.
1. Singularity (version 3.6.4-1.el7 or higher)
2. Docker (version 20.10.0 or higher)

[Download singularity here](https://singularity.lbl.gov/)  
[Download docker here](https://www.docker.com/products/docker-desktop)

Git and Conda are required to download the pipeline and create a new Conda environment with the required software

1. Start by cloning the GitHub repo  
    a. `git clone https://github.com/JoshLoecker/pipeline`
   b. This will create a new folder `pipeline` in your current directory. Move this to the folder you would like it to be contained in
2. Within the `pipeline` folder, navigate to the `setup` folder.
3. Run the `setup.sh` script by calling `./setup.sh`
