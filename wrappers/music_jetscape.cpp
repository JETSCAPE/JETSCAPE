// Copyright @ Bjoern Schenke, Sangyong Jeon, Charles Gale, and Chun Shen
#include <stdio.h>
#include <sys/stat.h>

#include <cstring>
#include "../src/fluid_dynamics.h"
#include "./music_jetscape.h"

using namespace std;


MPI_MUSIC::MPI_MUSIC() {
    hydro_status = NOT_START;
}


MPI_MUSIC::~MPI_MUSIC() {
    if (hydro_status != NOT_START) {
        delete music_hydro_ptr;
    }
}

void MPI_MUSIC::initialize_hydro(Parameter parameter_list) {
    int argc = 2;
    char **argv;
    music_hydro_ptr = new MUSIC(argc, argv);
    music_hydro_ptr->initialize_hydro();
    hydro_status = INITIALIZED;
}


void MPI_MUSIC::evolve_hydro() {
    music_hydro_ptr->run_hydro();
}



