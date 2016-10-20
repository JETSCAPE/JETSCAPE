// Copyright JETSCAPE Collaboration @ 2016

// This is a test script for JETSCAPE framework code to link with a 
// hydrodynamic code MUSIC
// Written by Chun Shen

#include <iostream>
#include <cstring>

#include "../src/fluid_dynamics.h"
#include "../wrappers/hydro_file_jetscape.h"

using namespace std;

int main(int argc, char *argv[]) {
    Parameter parameter_list;
    parameter_list.hydro_input_filename = *(argv+1);

    HydroFile *hydro_from_file_ptr = new HydroFile(parameter_list);
    hydro_from_file_ptr->initialize_hydro(parameter_list);
    hydro_from_file_ptr->evolve_hydro();

    return(1);
}
