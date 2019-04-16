#clone the urqmd-afterburner package from Jonah Bernhard
#git clone https://github.com/jbernhard/urqmd-afterburner.git

#clone version with modifications to osc2u from Chun (iEBE)
git clone https://github.com/derekeverett/urqmd-afterburner.git

#compile it
cd urqmd-afterburner

mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$WORK/software/ ..
make -j20 install

cd ..
