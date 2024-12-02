#!/usr/bin/env bash

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

# download the code package
#git clone git://git.code.sf.net/p/music-hydro/code music
#git clone --depth=1 git://git.code.sf.net/p/music-hydro/code music
#git clone --depth=1 https://github.com/MUSIC-fluid/MUSIC -b JETSCAPE music

# using a commit from the MUSIC repository that is compatible with JETSCAPE 3.6.6 and 4.0
folderName="music"
commitHash="fdcce5d9f63fc25281eda8d498c32c69ac132dac"

git clone https://github.com/MUSIC-fluid/MUSIC -b JETSCAPE $folderName
cd $folderName
git checkout $commitHash

### ALTERNATIVE VERSION
### Download a zipped snapshot
### This seems to require somebody visiting the download link on an actual browser
### to create a then-cached snapshot. That in turn seems to disappear about every day
### So revert back to using git (but make a shallow clone)

# # Not sure if we want to delete any MUSIC version already there...
# rm -rf music
# # solution from JF
# curlcmd1=wget
# curlcmd2="wget -qO-"
# command -v ${curlcmd1} > /dev/null || { curlcmd1="curl -LOk"; curlcmd2="curl -sk"; }
# command -v ${curlcmd1} > /dev/null || { echo "Please install curl or wget" ; exit 1; }

# downloadpath=$($curlcmd2 https://sourceforge.net/p/music-hydro/code/ci/public_release_prep/tarball | grep -o "https://sourceforge.net/[a-zA-Z0-9_/.-]*/music-hydro-code-[a-zA-Z0-9]*.zip" | head -1)
# filename="${downloadpath##*/}"
# filenameWithoutExt="${filename%%.*}"

# # Download file with either curl or wget
# $curlcmd1 ${downloadpath}

# # Put file in right location
# unzip ${filename} > /dev/null
# rm -f ${filename}
# mv ${filenameWithoutExt} music
