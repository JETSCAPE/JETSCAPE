#!/usr/bin/env bash

# download the code package
git clone https://git.code.sf.net/p/music-hydro/code music

# setup enviroment variables
echo "setup enviroment variables for MUSIC ..."
export HYDROPROGRAMPATH="${PWD}/music"
echo "HYDROPROGRAMPATH = $HYDROPROGRAMPATH"
