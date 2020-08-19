# Pipeline
## Installation and Running
This project is meant to run using Singularity or Docker
### <ins>Folder Setup</ins>
Two paths are needed with this container: 1) A `results` folder, and 2) A `data_files` folder. These folders can be named as you please, but this guide will use these respective names
It is important to note that under the `data_files` folder, a folder **must** be named `fast5`. If this is not done, guppy_basecaller will not be able to pick up **any** `.fast5` files, even if they are present, as it looks for `.fast5` files under the `fast5` folder. As such, your `data_files` folder structure may look as follows:
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
The `results` folder does not need to exist, but it can if you would like. If it does not exist, It will be created by Singularity/Docker when starting the container
### <ins>Singularty</ins>
1. If you are working on a cluster, Singularity may already be installed; this is the case with SciNet. If Singularity is not installed, follow the instructions: [Install Singularity](https://singularity.lbl.gov/install-linux)
2. Pull the container by executing the following command in a terminal window:<br>
	`singularity pull docker://joshloecker/pipeline:latest`
	
	a) This will begin downloading the Docker container and converting it to a Singularity container
4. To run the container, execute the following command:
    ```
    singularity run \
    -B /path/to/results/:/results/ \
	-B /path/to/data_files/:/data_files/ \
	pipeline:latest
	```
### <ins>Docker</ins>
1. If you chose to work with Docker, ensure Docker is already installed on your system. If it is not, follow the instructions: [Install Docker](https://docs.docker.com/get-docker/)
2. Pull the container by executing the following command in a terminal window:<br>
	`docker pull joshloecker/pipeline:latest`
3. Run the container using the following command:
	```
	docker run \
	-v /path/to/results/:/results/ \
	-v /path/to/data_files/:/data_files/ \
 	pipeline:latest
 	```

 ### <ins>Known Errors</ins>
 This pipeline is not 100% complete, as the ability to bind an alignment file needs to be added.

 If you experience any errors, contact joshua.loecker@usda.gov for assistance


 ## Results
 Results from the pipeline will be saved to the path you bound during the **Installation and Running** step (i.e. `/path/to/results/`). From here, you will be able to see each of the folders that are output, such as `basecalling`, `barcoding`, `visuals`, etc.
