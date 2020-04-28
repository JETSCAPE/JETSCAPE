# Build from the official docker python base image, based on Debian
FROM python:3.6-stretch

# Install pre-reqs (commented ones are already in base image)
RUN apt-get update && apt-get install -y \
cmake \
#curl \
doxygen \
emacs \
#g++ \
gsl-bin \
hdf5-tools \
less \
libboost-all-dev \
libeigen3-dev \
libgsl-dev \
libhdf5-serial-dev \
libxpm-dev \
openmpi-bin \
rsync \
valgrind \
vim \
#zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*

# Install additional useful python packages
RUN pip install \
pyyaml

# Install ROOT6 from source
ENV ROOTSYS="/usr/local/root"
ENV PATH="${ROOTSYS}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${ROOTSYS}/lib:${LD_LIBRARY_PATH}"
ENV PYTHONPATH="${ROOTSYS}/lib"
RUN mkdir -p ${ROOTSYS} && mkdir -p ${HOME}/root && cd ${HOME}/root \
&& git clone --depth 1 --branch v6-14-06 https://github.com/root-project/root.git src \
&& mkdir build && cd build \
&& cmake ../src -DCMAKE_INSTALL_PREFIX=${ROOTSYS} \
&& make -j8 install \
&& rm -r ${HOME}/root

# Install HepMC 3.1.1
RUN curl -SL http://hepmc.web.cern.ch/hepmc/releases/HepMC3-3.1.1.tar.gz | tar -xvzC /usr/local \
&& cd /usr/local \
&& mkdir hepmc3-build \
&& cd hepmc3-build \
&& cmake ../HepMC3-3.1.1 -DCMAKE_INSTALL_PREFIX=/usr/local -DHEPMC3_ENABLE_ROOTIO=ON -DROOT_DIR=${ROOTSYS} -DHEPMC3_BUILD_EXAMPLES=ON \
&& make -j8 install \
&& rm -r /usr/local/hepmc3-build
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

# Install Pythia8 (8.235 is needed by SMASH, and it should be installed in a folder called pythia8235)
RUN curl -SL http://home.thep.lu.se/~torbjorn/pythia8/pythia8235.tgz \
| tar -xvzC /usr/local \
&& cd /usr/local/pythia8235 \
&& ./configure --enable-shared --prefix=/usr/local/pythia8235 --with-hepmc3=/usr/local/HepMC3-3.1.1 \
&& make -j8 \
&& make -j8 install

# Install FastJet
RUN curl -SL http://fastjet.fr/repo/fastjet-3.3.2.tar.gz \
| tar -xvzC /usr/local \
&& cd /usr/local/fastjet-3.3.2 \
&& ./configure --prefix=/usr/local \
&& make -j8 \
&& make -j8 install \
&& rm -r /usr/local/fastjet-3.3.2

# Install FastJet Contrib
RUN curl -SL http://fastjet.hepforge.org/contrib/downloads/fjcontrib-1.041.tar.gz \
| tar -xvzC /usr/local \
&& cd /usr/local/fjcontrib-1.041 \
&& ./configure \
&& make -j8 \
&& make -j8 install \
&& rm -r /usr/local/fjcontrib-1.041

# Set environmental variables for cmake to know where things are (needed for SMASH)
ARG username=jetscape-user
ENV JETSCAPE_DIR="/home/${username}/JETSCAPE"
ENV SMASH_DIR="${JETSCAPE_DIR}/external_packages/smash/smash_code"
ENV EIGEN3_ROOT /usr/include/eigen3
ENV PYTHIA8DIR /usr/local/pythia8235
ENV PYTHIA8_ROOT_DIR /usr/local/pythia8235

# Set up a user group
ARG id=1234
RUN groupadd -g ${id} ${username} \
&& useradd --create-home --home-dir /home/${username} -u ${id} -g ${username} ${username}
USER ${username}
ENV HOME /home/${username}
WORKDIR ${HOME}

ENTRYPOINT /bin/bash
