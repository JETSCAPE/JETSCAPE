#!/usr/bin/env bash

# download the code package
#git clone git://git.code.sf.net/p/music-hydro/code music

# solution from JF
curlcmd=wget
command -v ${curlcmd} > /dev/null || curlcmd="curl -LO"
command -v ${curlcmd} > /dev/null || { echo "Please install curl or wget" ; exit 1; }
wget $(wget -qO- https://sourceforge.net/p/music-hydro/code/ci/public_release_prep/tarball | grep -o "https://sourceforge.net/[a-zA-Z0-9_/.-]*/music-hydro-code-[a-zA-Z0-9]*.zip" | head -1)
unzip music-hydro-code-[a-zA-Z0-9]*.zip > /dev/null
rm -fr music-hydro-code-[a-zA-Z0-9]*.zip
ii=`ls | grep "music-hydro-code-[a-zA-Z0-9]*"`
mv $ii music
