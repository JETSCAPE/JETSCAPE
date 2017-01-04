// Copyright JETSCAPE Collaboration @ 2016

// This is a test script for JETSCAPE framework code to link with a 
// hydrodynamic code MUSIC
// Written by Chun Shen

#include <iostream>
#include <cstring>

#include "fluid_dynamics.h"
#include "hydro_file_jetscape.h"

using namespace std;

int main(int argc, char *argv[]) {
    Parameter parameter_list;
    parameter_list.hydro_input_filename = *(argv+1);

    HydroFile *hydro_from_file_ptr = new HydroFile();
    hydro_from_file_ptr->initialize_hydro(parameter_list);
    hydro_from_file_ptr->evolve_hydro();
    FluidCellInfo* check_fluid_info_ptr = new FluidCellInfo;
    hydro_from_file_ptr->get_hydro_info(1.0, 0.0, 0.0, 0.0,
                                        check_fluid_info_ptr);
    hydro_from_file_ptr->print_fluid_cell_information(check_fluid_info_ptr);
    return(1);
}
