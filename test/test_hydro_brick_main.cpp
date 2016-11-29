// Copyright JETSCAPE Collaboration @ 2016

// This is a test script for JETSCAPE framework code to link with hydro brick
// Written by Chun Shen

#include <iostream>
#include <cstring>

#include "../src/fluid_dynamics.h"
#include "../wrappers/hydro_brick_jetscape.h"

using namespace std;

int main(int argc, char *argv[]) {
    Parameter parameter_list;
    parameter_list.hydro_input_filename = *(argv+1);

    HydroBrick *hydro_brick_ptr = new HydroBrick();
    hydro_brick_ptr->initialize_hydro(parameter_list);
    hydro_brick_ptr->evolve_hydro();
    FluidCellInfo* check_fluid_info_ptr = new FluidCellInfo;
    hydro_brick_ptr->get_hydro_info(1.0, 0.0, 0.0, 0.0, check_fluid_info_ptr);
    hydro_brick_ptr->print_fluid_cell_information(check_fluid_info_ptr);
    return(0);
}
