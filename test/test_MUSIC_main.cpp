#include <iostream>
#include <cstring>

#include "../src/fluid_dynamics.h"
#include "./music.h"

using namespace std;

int main(int argc, char *argv[]) {
    Parameter parameter_list;
    parameter_list.hydro_input_filename = *(argv+1);
    MUSIC *MUSIC_ptr = new MUSIC();
    MUSIC_ptr->initialize_hydro(parameter_list);
    MUSIC_ptr->evolve_hydro();
}
