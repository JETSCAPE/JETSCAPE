***The JETSCAPE / SIMS branch*** : 

Currently, the SIMS branch integrates up-to-date bulk medium simulation modules:

1. Trento initial condition (2d and 3d)
2. 3+1D freestreaming
3. 3+1D relativitic viscous hydrodynamics MUSIC
4. iS3D particlization module
5. The hadronic transport model SMASH

# Load libraries on stampede2

   ```bash
      module load intel/18.0.2 cmake/3.7.1 gsl boost hdf5 eigen impi
   ```

# Install Pythia8235

   ```bash
      wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8235.tgz
      tar xvfz pythia8235.tgz && rm pythia8235.tgz
      cd pythia8235
   ```

   Configure it to install into $WORK/software

   ```bash
      ./configure --prefix=$WORK/software
      make
      make install
   ```

   Setup enviroment for pythia

   ```bash
      echo "export PATH=\$WORK/software/bin:\$PATH" >> $HOME/.bashrc
      echo "export LD_LIBRARY_PATH=\$WORK/software/lib:\$LD_LIBRARY_PATH" >> $HOME/.bashrc
      echo "export PYTHIA8DATA=\$WORK/software/share/Pythia8/xmldoc/" >> $HOME/.bashrc
      echo "export PYTHIA8=\$HOME/software" >> $HOME/.bashrc
      source $HOME/.bashrc
   ```


# Compilation on stampede2 TACC

1. Clone the SIMS branch into the $WORK directory (make sure it is sims branch)

   ```bash
      cd $WORK && pwd
      git clone -b sims https://github.com/amajumder/JETSCAPE-COMP.git
      cd JETSCAPE-COMP && git branch
   ```

2. Remark: SMASH has to be compiled before JETSCAPE. 
  
   Stempede has cmake, gsl, boost and eigen3 installed already, except for pythia which should have been in you system path from the last section. Then, activate all the environments by running the sims activation script (make sure the paths are set correctly)

   ```bash
      source <JETSCAPE>/sims-scripts/prepare_compilation_stampede2_2.sh
   ```

   Then, go to the external modules directory and download smash (also compiles smash), iS3D, freestream-milne, music by
   
   ```bash
      cd <JETSCAPE>/external_packages/
      ./get_modules_for_sims.sh
   ```

3. Compile JETSCAPE SIMS program

   Go back to the JETSCAPE folder, and then

   ```bash
      mkdir build && cd build
      cmake -Dmusic=on -DiS3D=on -Dsmash=on -Dfreestream=on ../
   ```
   
   In case you want to install the executables to a designated path and later add to system path (e.g. $WORK/software/bin). You can add the path flag to the Cmake command by,
   
   ```bash
      cmake -DCMAKE_INSTALL_PREFIX=$WORK/software/ -Dmusic=on -DiS3D=on -Dsmash=on -Dfreestream=on ../ 
   ```

   Remember to switch on all those compilation flags to use external modules.

   ```bash
      make
      make install
   ``` 

4. Run simulation

   Before running the program somewhere else than the build folder, or in the submit scrit. Makesure you have set the followings

   ```bash
      export PATH=$WORK/software/bin:$PATH
      export LD_LIBRARY_PATH=$WORK/software/lib:$LD_LIBRARY_PATH
      export PYTHIA8DATA=$WORK/software/share/Pythia8/xmldoc
      export PYTHIA8=$HOME/software
   ```

   Also remember to copy the inputfiles (can be found in build/ after compilation) for all programs to that folder, including

   * `EOS/` : contains hydro equation of state tables
   * `freestream_input` : freestream configure file
   * `iS3D_parameters.dat` : particle sampler configure file
   * `music_input` : hydro configure file
   * `smash_config/` : a folder containing smash config data
   * `jetscape_init.xml` : the JetScape configure file.

   Remark: initial condition normalization should be set by `jetscape_init.xml` via the `trento` block, the `s_factor` in `music_input` should be set to `1.0`.

   Now it should be good to go. The first executable `TRENTO_FS_HYDRO` runs trento initial condition, freestream and music hydro. The second executable `SAMPLER_AFTERBURNER` runs particle sampler iS3D and SMASH hadronic afterburner. 

5. Todo
 
   * A workflow control python script is in progress
   * Provide sample configure files
   * Scripts for calculating observables for each event
   * Sample job submisstion script on stampede2
   * and more...

