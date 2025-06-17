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

# This script runs a Docker container for a specified version of
# JETSCAPE/X-SCAPE.  Docker must be installed on the host system
# and accessible from the Linux or WSL bash shell.

# This script takes two command line arguments:
# 1) The path to the user input XML file.
# 2) The image tag for the JETSCAPE version, e.g. beta_v0.11

# The available image tags for JETSCAPE are listed at:
# https://hub.docker.com/r/jetscape/jetscape_full

# Example usage:
# ./runDocker.sh arXiv_1910.05481/jetscape_user_PP_1910.05481.xml beta_v0.11

# check command line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_xml> <image_repository:tag>"
    echo "Example: $0 jetscape_user.xml jetscape_full:beta_v0.11"
    exit 1
fi

# set accordingly for JETSCAPE or X-SCAPE
repo_name="JETSCAPE"

# set accordingly for jetscape_full or xscape_full
image_repo_name="jetscape_full"

input_xml=$1
image_tag=$2

# check if input XML file exists
if [ ! -f "$input_xml" ]; then
    echo "Error: Input XML file '$input_xml' not found!"
    exit 1
fi

build_dir="/home/jetscape-user/$repo_name/build"
entry_point="/home/jetscape-user/$repo_name/build/runJetscape"
host_mount_dir="/home/jetscape-user/$repo_name/host"

# run the Docker container
docker run --rm \
    -w "$build_dir" \
    --user "$(id -u):$(id -g)" \
    --entrypoint "$entry_point" \
    -v $(pwd):"$host_mount_dir" \
    "jetscape/$image_repo_name:$image_tag" \
    "$host_mount_dir/$input_xml"
