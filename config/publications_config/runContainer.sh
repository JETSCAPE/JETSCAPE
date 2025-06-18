#!/bin/bash

###############################################################################
# Copyright (c) The JETSCAPE Collaboration, 2018
#
# For the list of contributors see AUTHORS.
#
# Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
#
# or via email to bugs.jetscape@gmail.com
#
# Distributed under the GNU General Public License 3.0 (GPLv3 or later).
# See COPYING for details.
##############################################################################

# This script runs a Docker or Apptainer container for a specified
# version of JETSCAPE or X-SCAPE. Docker or Apptainer must be installed on
# the host system and accessible from the Linux bash shell.

# This script takes two command line arguments:
# 1) The path to the user input XML file.
# 2) The image repository and tag for the JETSCAPE or X-SCAPE version.
#    For example: jetscape_full:v3.7.1

# The JETSCAPE and X-SCAPE images are available at:
# https://hub.docker.com/r/jetscape/jetscape_full
# https://hub.docker.com/r/jetscape/xscape_full

# Example usage:
# ./runContainer.sh arXiv_1910.05481/jetscape_user_PP_1910.05481.xml jetscape_full:v3.7.1

# check command line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_xml> <image_repository:tag>"
    echo "Example: $0 arXiv_1910.05481/jetscape_user_PP_1910.05481.xml jetscape_full:v3.7.1"
    exit 1
fi

input_xml=$1
image_repo_tag=$2

# set accordingly for JETSCAPE or X-SCAPE
if [[ "$image_repo_tag" == jetscape_full* ]]; then
    repo_name="JETSCAPE"
elif [[ "$image_repo_tag" == xscape_full* ]]; then
    repo_name="X-SCAPE"
else
    echo "Error: image repository name must be jetscape_full or xscape_full"
    exit 1
fi

# check if input XML file exists
if [ ! -f "$input_xml" ]; then
    echo "Error: Input XML file '$input_xml' not found"
    exit 1
fi

build_dir="/home/jetscape-user/$repo_name/build"
entry_point="/home/jetscape-user/$repo_name/build/runJetscape"
host_mount_dir="/home/jetscape-user/$repo_name/host"

# run the Docker container if Docker is available
if command -v docker &> /dev/null; then
    echo "Running Docker command..."
    docker run --rm \
        -w "$build_dir" \
        --user "$(id -u):$(id -g)" \
        --entrypoint "$entry_point" \
        -v $(pwd):"$host_mount_dir" \
        "jetscape/$image_repo_tag" \
        "$host_mount_dir/$input_xml"

# run the Apptainer container if Apptainer is available
elif command -v apptainer &> /dev/null; then
    echo "Running Apptainer command..."
    apptainer exec \
        --pwd "$build_dir" \
        --bind "$(pwd):$host_mount_dir" \
        docker://jetscape/$image_repo_tag \
        "$entry_point" "$host_mount_dir/$input_xml"

# run the Singularity container if Singularity is available
elif command -v singularity &> /dev/null; then
    echo "Running Singularity command..."
    singularity exec \
        --pwd "$build_dir" \
        --bind "$(pwd):$host_mount_dir" \
        docker://jetscape/$image_repo_tag \
        "$entry_point" "$host_mount_dir/$input_xml"

# error if none are available
else
    echo "Error: neither Docker nor Apptainer is available"
    exit 1
fi
