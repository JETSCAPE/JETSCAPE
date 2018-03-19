#!/usr/bin/env bash

# download the code package
#git clone git://git.code.sf.net/p/music-hydro/code music

# Not sure if we want to delete any MUSIC version already there...
rm -rf music
# solution from JF
curlcmd1=wget
curlcmd2="wget -qO-"
command -v ${curlcmd1} > /dev/null || { curlcmd1="curl -LOk"; curlcmd2="curl -sk"; }
command -v ${curlcmd1} > /dev/null || { echo "Please install curl or wget" ; exit 1; }

downloadpath=$($curlcmd2 https://sourceforge.net/p/music-hydro/code/ci/public_release_prep/tarball | grep -o "https://sourceforge.net/[a-zA-Z0-9_/.-]*/music-hydro-code-[a-zA-Z0-9]*.zip" | head -1)
filename="${downloadpath##*/}"
filenameWithoutExt="${filename%%.*}"

# Download file with either curl or wget
$curlcmd1 ${downloadpath}

# Put file in right location
unzip ${filename} > /dev/null
rm -f ${filename}
mv ${filenameWithoutExt} music
