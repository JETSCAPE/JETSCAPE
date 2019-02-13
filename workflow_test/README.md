# ***The JETSCAPE / SIMS branch*** : 

Currently, the SIMS branch integrates up-to-date bulk medium simulation modules:

1. Trento initial condition (2d and 3d)
2. 3+1D freestreaming
3. 3+1D relativitic viscous hydrodynamics MUSIC
4. iS3D particlization module
5. The hadronic transport model SMASH

# Quick Install

The fastest way to install the sims branch on stampede2 is to use the workflow test branch. Downloading `sims_install.sh` and copying it to stampede2 is the cleanest way.
First:
   ```bash
      sh sims_install.sh
   ```
This requires GitHub credentials, but otherwise can be left alone. The process takes approx. 15-20 minutes.

The most straightforward way is to:
   ```bash
      cd $WORK
      git clone -b sims https://github.com/amajumder/JETSCAPE-COMP.git 
      cd workflow_test/
   ```
***Then comment out line 12 which clones the Git repository***

Finally,
   ```bash
      sh sims_install.sh
   ```

To test the installation and workflow,
   ```bash
      sh sims_workflow.sh
   ```
 ***Remark***: initial condition normalization parameter should be set in `jetscape_init.xml` via the `trento` block, the `s_factor` in `music_input` should be set to `1.0`. If `sims_workflow.sh` encounters difficulty, check the `s_factor`.

This process follows the instructions in the sims scripts readme, but automates the process in detail and is self-contained, so is unlikely to be affected by changes to the latest version of the workflow.

However, if modifying components and recompiling outside the workflow test script afterward, several of the steps below will be necessary. 

1. Load libraries on stampede2
2. Setup environment for pythia (at least ensure paths are exported to .bashrc)

***Remark***: The following must be done before recompiling and running each time
3. Proceed from step 2 of "Compilation on stampede2 TACC" below.

Additional steps that must be taken can be found by following `sims_workflow.sh`, but explicitly:
1. Make a directory for running an event, say you call it `test_newworkflow`
2. From `workflow_test`, copy `generate_module_input_files.py`, `Run_WorkflowTest.py`, `workflow_test`, `workflowtest_oversampled_afterburner.py`, `smash_input/` into `test_newworkflow`
3. From `/external_packages/iS3D/` copy  `PDG/`, `tables/`, `deltaf_coefficients/` into `test_newworkflow`
4. From `/external_packages/music/` copy ` EOS/` into your `test_worflow` directory
Then `python generate_module_input_files.py`  and   `sbatch workflow_test`

# Full install instructions

## Load libraries on stampede2

   ```bash
      module load intel/18.0.2 cmake/3.7.1 gsl boost hdf5 eigen impi
   ```

## Install Pythia8235

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
      echo "export PYTHIA8=\$WORK/software" >> $HOME/.bashrc
      source $HOME/.bashrc
   ```


## Compilation on stampede2 TACC

1. Clone the SIMS branch into the `$WORK` directory (make sure it is sims branch)

   ```bash
      cd $WORK && pwd
      git clone -b sims https://github.com/amajumder/JETSCAPE-COMP.git
      cd JETSCAPE-COMP && git branch
   ```

2. Remark: SMASH has to be compiled before JETSCAPE. 
  
   Stempede has cmake, gsl, boost and eigen3 installed already, except for pythia which should have been in you system path from the last section. Then, activate all the environments by running the sims activation script (make sure the paths are set correctly)

   ```bash
      source <JETSCAPE>/sims_scripts/prepare_compilation_stampede2_2.sh
   ```

   Then, go to the external modules directory and download `smash` (also compiles `smash`), `iS3D`, `freestream-milne`, `music` by
   
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

## Run simulation

   If you want to run the program elsewhere than the build folder, or if you are writting the job submit scrit. Make sure you have set the `PATH`, `LD_LIBRARY_PATH` and PYTHIA paths. For example, create a testing folder under `build/`
   
   ```bash
      mkdir <build>/test
      cd <build>/test
   ```
   
   Then, set the following to have things work properly, 

   ```bash
      export PATH=$WORK/software/bin:$PATH
      export LD_LIBRARY_PATH=$WORK/software/lib:$LD_LIBRARY_PATH
      export PYTHIA8DATA=$WORK/software/share/Pythia8/xmldoc
      export PYTHIA8=$WORK/software
   ```

   Also, remember to copy the inputfiles (can be found in `build/` after compilation) for all programs to the test folder, including,
   
   * `freestream_input` : free-stream configure file
   * `EOS/` : contains hydro equation-of-state tables
   * `music_input` : hydro configure file  
   * `iS3D_parameters.dat` : particle sampler iS3D configure file
   * `deltaf_coefficients/` : contains delta-f correction coefficients
   * `PDG/` : hadron info for sampler iS3D
   * `tables` : integration tables for sampler iS3D
   * `smash_config/` : a folder that contains smash configure files
   * `jetscape_init.xml` : the JETSCAPE configure file

   ***Remark***: initial condition normalization parameter should be set in `jetscape_init.xml` via the `trento` block, the `s_factor` in `music_input` should be set to `1.0`.

   * Also ensure that `generate_module_input_files.py` has been copied and `python generate_module_input_files.py`

   Now it should be good to go. The first executable `TRENTO_FS_HYDRO` runs trento initial condition, freestream and music hydro. Then `mkdir input && cp surface.dat input/` before running the next step. The second executable `SAMPLER_AFTERBURNER` runs particle sampler iS3D and SMASH hadronic afterburner. 

## Todo
 
      * We should provide a set of sample configure files
   * A script to calculate observables for each event
   * A sample job submisstion script for stampede2
   * and more...

