// Copyright @ Bjoern Schenke, Sangyong Jeon, Charles Gale, and Chun Shen
#include <stdio.h>
#include <sys/stat.h>

#include <cstring>

#include "JetScapeLogger.h"
#include "music_jetscape.h"

using namespace std;

MPI_MUSIC::MPI_MUSIC() {
    hydro_status = NOT_START;
    SetId("MUSIC");
}


MPI_MUSIC::~MPI_MUSIC() {
    if (hydro_status != NOT_START) {
        delete music_hydro_ptr;
    }
}


void MPI_MUSIC::initialize_hydro(Parameter parameter_list) {
    INFO << "Initialize MUSIC ...";
    VERBOSE(8);
    tinyxml2::XMLElement *para =
                    GetHydroXML()->FirstChildElement("MUSIC");
    if (!para) {
        WARN << " : MUSIC not properly initialized in XML file ...";
        exit(-1);
    }
    string input_file = para->FirstChildElement("MUSIC_input_file")->GetText();
    int argc = 2;
    char **argv = new char* [argc];
    argv[0] = new char[8];
    strcpy(argv[0], "mpihydro");
    argv[1] = new char[input_file.length() + 1];
    strcpy(argv[1], input_file.c_str());
    cout << "check input for MUSIC: " << endl;
    for (int i = 0; i < argc; i++) {
        cout << argv[i] << "  ";
    }
    cout << endl;
    music_hydro_ptr = new MUSIC(argc, argv);
    //music_hydro_ptr->initialize_hydro();
    //hydro_status = INITIALIZED;

    for (int i = 0; i < argc; i++) {
        delete[] argv[i];
    }
    delete[] argv;
}


void MPI_MUSIC::evolve_hydro() {
    VERBOSE(8);
    std::vector<double> entropy_density = ini->entropy_density_distribution_;
    music_hydro_ptr->initialize_hydro_from_vector(entropy_density);
    hydro_status = INITIALIZED;
    if (hydro_status == INITIALIZED) {
        INFO << "running MUSIC ...";
        music_hydro_ptr->run_hydro();
        hydro_status = FINISHED;
    } else if (hydro_status == FINISHED) {
        INFO << "Initialize MUSIC ...";
        music_hydro_ptr->initialize_hydro();
        INFO << "running MUSIC ...";
        music_hydro_ptr->run_hydro();
        hydro_status = FINISHED;
    }
}


void MPI_MUSIC::get_hydro_info(real t, real x, real y, real z,
                               FluidCellInfo* fluid_cell_info_ptr) {
    fluidCell *fluidCell_ptr = new fluidCell;
    music_hydro_ptr->get_hydro_info(x, y, z, t, fluidCell_ptr);
    fluid_cell_info_ptr->energy_density = fluidCell_ptr->ed;
    fluid_cell_info_ptr->entropy_density = fluidCell_ptr->sd;
    fluid_cell_info_ptr->temperature = fluidCell_ptr->temperature;
    fluid_cell_info_ptr->pressure = fluidCell_ptr->pressure;
    fluid_cell_info_ptr->vx = fluidCell_ptr->vx;
    fluid_cell_info_ptr->vy = fluidCell_ptr->vy;
    fluid_cell_info_ptr->vz = fluidCell_ptr->vz;
    fluid_cell_info_ptr->mu_B = 0.0;
    fluid_cell_info_ptr->mu_C = 0.0;
    fluid_cell_info_ptr->mu_S = 0.0;
    fluid_cell_info_ptr->qgp_fraction = 0.0;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            fluid_cell_info_ptr->pi[i][j] = fluidCell_ptr->pi[i][j];
        }
    }
    fluid_cell_info_ptr->bulk_Pi = fluidCell_ptr->bulkPi;

    delete fluidCell_ptr;
}

