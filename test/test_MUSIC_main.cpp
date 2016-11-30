// Copyright JETSCAPE Collaboration @ 2016

// This is a test script for JETSCAPE framework code to link with a 
// hydrodynamic code MUSIC
// Written by Chun Shen

#include <iostream>
#include <cstring>

#include "../src/fluid_dynamics.h"
#include "../wrappers/music_jetscape.h"

using namespace std;

int main(int argc, char *argv[]) {
    Parameter parameter_list;
    parameter_list.hydro_input_filename = *(argv+1);

    MPI_MUSIC *MUSIC_ptr = new MPI_MUSIC();
    MUSIC_ptr->initialize_hydro(parameter_list);
    MUSIC_ptr->evolve_hydro();

    return(1);
}
