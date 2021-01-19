#!/usr/bin/env bash
# This script will perform the required steps to build a new Guppy container for an updated image
# Before running this script, be sure to update the `version` argument in the Dockerfile

# Set the guppy version you would like to use
# Visit https://community.nanoporetech.com/downloads to see what the latest version is
GUPPY_VERSION=4.4.1

# build guppy container
printf "Building guppy container\n"
docker build --build-arg version=$GUPPY_VERSION --tag guppy_container:latest .

docker_exit_code=$?
if [[ $docker_exit_code -gt 0 ]]; then
    echo "Start docker before continuing"
    echo "EXIT: $docker_exit_code"
    exit $docker_exit_code
fi

# build singularity image
printf "Building singularity container\n"
singularity build --sandbox guppy_container docker-daemon://guppy_container:latest
singularity_exit_code=$?
if [[ $singularity_exit_code -gt 0 ]]; then
    echo "Start singularity before continuing"
    echo "EXIT: $singularity_exit_code"
    exit $singularity_exit_code
fi

printf "Pruning docker system\n"
# prune docker of un-needed images, containers, etc.
docker system prune
