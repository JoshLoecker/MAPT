#!/usr/bin/env bash
# This script will perform the required steps to build a new Guppy container for an updated image
# Before running this script, be sure to update the `version` argument in the Dockerfile

# go to parent directory (container Dockerfile)
cd ../

# build guppy container
printf "Building guppy container\n"
docker build --tag guppy_container:latest .

# build singularity image
printf "Building singularity container\n"
singularity build --sandbox guppy_container docker-daemon://guppy_container:latest

printf "Pruning docker system\n"
# prune docker of un-needed images, containers, etc.
docker system prune
