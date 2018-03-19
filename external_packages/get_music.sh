#!/usr/bin/env bash

# download the code package
#git clone git://git.code.sf.net/p/music-hydro/code music

rm -rf music
rm -rf music-hydro-code-[a-zA-Z0-9]*
# solution from JF
curlcmd1=wget
curlcmd2="wget -qO-"
command -v ${curlcmd1} > /dev/null || { curlcmd1="curl -LOk"; curlcmd2="curl -sk"; }
command -v ${curlcmd1} > /dev/null || { echo "Please install curl or wget" ; exit 1; }

$curlcmd1 $($curlcmd2 https://sourceforge.net/p/music-hydro/code/ci/public_release_prep/tarball | grep -o "https://sourceforge.net/[a-zA-Z0-9_/.-]*/music-hydro-code-[a-zA-Z0-9]*.zip" | head -1)

unzip music-hydro-code-[a-zA-Z0-9]*.zip > /dev/null
rm -fr music-hydro-code-[a-zA-Z0-9]*.zip
ii=`ls | grep "music-hydro-code-[a-zA-Z0-9]*"`
mv $ii music
