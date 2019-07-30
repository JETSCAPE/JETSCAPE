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
libgsl-dev \
libhdf5-serial-dev \
libxpm-dev \
openmpi-bin \
rsync \
vim \
#zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*

# Install a few useful python packages
RUN pip3 install \
pyyaml
#subprocess

# Install HEPMC3
RUN curl -SL http://hepmc.web.cern.ch/hepmc/releases/hepmc3.0.0.tgz \
| tar -xvzC /usr/local \
&& cd /usr/local/hepmc3.0.0/cmake \
&& cmake .. \
&& make -j8 all install \
&& rm -r /usr/local/hepmc3.0.0

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

# Install Pythia8
RUN curl -SL http://home.thep.lu.se/~torbjorn/pythia8/pythia8240.tgz \
| tar -xvzC /usr/local \
&& cd /usr/local/pythia8240 \
&& ./configure --enable-shared --prefix=/usr/local --with-hepmc3=/user/local/hepmc3.0.0 \
&& make -j8 \
&& make -j8 install \
&& rm -r /usr/local/pythia8240

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
# Note: The build is not cached.

# Set up a user group
ARG username=jetscape-user
ARG id=1234
RUN groupadd -g ${id} ${username} \
&& useradd --create-home --home-dir /home/${username} -u ${id} -g ${username} ${username}
USER ${username}
ENV HOME /home/${username}
WORKDIR ${HOME}

ENTRYPOINT /bin/bash
