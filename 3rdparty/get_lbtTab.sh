#!/usr/bin/env bash

# download the tables required by LBT code
curlcmd=wget
command -v ${curlcmd} > /dev/null || curlcmd="curl -LO"
command -v ${curlcmd} > /dev/null || { echo "Please install curl or wget" ; exit 1; }
$curlcmd "https://bitbucket.org/sscao/lbt-tables/downloads/LBT-tables.tar.gz"

tar xvzf LBT-tables.tar.gz

