#!/usr/bin/env bash
# This file will create a conda environment using the environment file in the parent directory
# If the conda environment already exists, nothing will be done

# Get the prefix name
# This is also the name of the environment
prefix_name=$(tail < "../environment.yaml" | grep "prefix: ")
prefix_name="${prefix_name:8}"  # take data after 'prefix: ' (i.e. prefix name)

# Get a list of current environments
current_conda_envs=$(conda env list)

# we have not found an environment with the name 'mapt_pipeline', we need to create it
if [[ "$current_conda_envs" != *"$prefix_name"* ]]; then
    # Create a new conda environment from the "environment.yaml" file in the parent directory
    conda env create --file "../environment.yaml" --prefix "${prefix_name}"
    # only use the environment name when activating, not the full prefix path
    conda config --set env_prompt '({name})'
# we have found a conda environment with the name 'pipeline', check if it has the correct packages
else
    printf "Found the required conda environment: %s\n" "$prefix_name"
fi

if [[ "$?" == 0 ]]; then
    printf "\n"
    printf "Setup complete.\n"
    printf "You can list all environments with 'conda info --envs'\n"
    printf "The MAPT pipeline will be using the environment %s\n" "$prefix_name"
    printf "Activate this environment with 'conda activate %s\n'" "$prefix_name"
fi
