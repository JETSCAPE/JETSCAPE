#!/usr/bin/env bash

# download VISHNU hydro examples
curlcmd=wget
command -v ${curlcmd} > /dev/null || curlcmd="curl -LO"
command -v ${curlcmd} > /dev/null || { echo "Please install curl or wget" ; exit 1; }
$curlcmd "https://bitbucket.org/sscao/hydro_profiles/downloads/PbPb2760-Avg-00-05.tar.gz"

tar xvzf PbPb2760-Avg-00-05.tar.gz
mkdir -p test_hydro_files
mv PbPb2760-Avg-00-05 test_hydro_files/event-0
rm PbPb2760-Avg-00-05.tar.gz
