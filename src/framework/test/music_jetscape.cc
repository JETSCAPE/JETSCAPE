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
    int argc = 2;
    char **argv = nullptr;
    // music_hydro_ptr = new MUSIC(argc, argv);
    // music_hydro_ptr->initialize_hydro();
    hydro_status = INITIALIZED;
}


void MPI_MUSIC::evolve_hydro() {
    VERBOSE(8);
    INFO << "running MUSIC ...";
    // music_hydro_ptr->run_hydro();
    hydro_status = FINISHED;
}



