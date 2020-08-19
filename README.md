# Pipeline

## Installation and Running
This project is meant to run using Singularity or Docker

### <ins>Singularty</ins>
1. If you are working on a cluster, Singularity may already be installed; this is the case with SciNet. If Singularity is not installed, follow the instructions: [Install Singularity](https://singularity.lbl.gov/install-linux)
2. Begin downloading the Pipeline Container by running `singularity pull docker://joshloecker/pipeline:latest`
	a) This will being downloading the Docker container and converting it to a Singularity container

Two paths are needed with this container: 1) A `results` folder, and 2) A `data_files` folder. These folders can be named as you please, but this guide will use these respective names

3. To run the container, execute the following command:

	```
	singularity run \
	-B /path/to/results/folder:/results/ \
	-B /path/to/data_files/folder:/data_files/ \
	pipeline:latest
	```


### <ins>Docker</ins>
1. If you chose to work with Docker, ensure Docker is already installed on your system. If it is not, follow the instructions: [Install Docker](https://docs.docker.com/get-docker/)



THIS IS SOME TEXT
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTc3MzI1OTA1N119
-->