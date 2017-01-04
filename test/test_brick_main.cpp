// Copyright JETSCAPE Collaboration @ 2016

// This is a test script for JETSCAPE framework code to link with brick
// Written by Chun Shen

#include <iostream>
#include <cstring>

#include "fluid_dynamics.h"
#include "brick_jetscape.h"

using namespace std;

int main(int argc, char *argv[]) {
    Parameter parameter_list;
    parameter_list.hydro_input_filename = *(argv+1);

    Brick *brick_ptr = new Brick();
    brick_ptr->initialize_hydro(parameter_list);
    brick_ptr->evolve_hydro();
    FluidCellInfo* check_fluid_info_ptr = new FluidCellInfo;
    brick_ptr->get_hydro_info(1.0, 0.0, 0.0, 0.0, check_fluid_info_ptr);
    brick_ptr->print_fluid_cell_information(check_fluid_info_ptr);
    return(0);
}
