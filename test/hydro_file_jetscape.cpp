// Original copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
// Massively cleaned up and improved by Chun Shen 2015-2016
#include <stdio.h>
#include <sys/stat.h>

#include <cstring>
#include <iostream>
#include "../src/fluid_dynamics.h"
#include "./hydro_file_jetscape.h"
#include "../external_packages/hydro_from_external_file/Hydroinfo_h5.h"
#include "../external_packages/hydro_from_external_file/Hydroinfo_MUSIC.h"
#include "../external_packages/hydro_from_external_file/ParameterReader.h"

using namespace std;


HydroFile::HydroFile(Parameter parameter_list) {
    // initialize the parameter reader
    paraRdr = new ParameterReader();
    paraRdr->readFromFile("parameters.dat");
    paraRdr->echo();

    hydro_status = 0;
}


HydroFile::~HydroFile() {
    if (hydro_status > 0) {
        delete paraRdr;
    }
}


void HydroFile::initialize_hydro(Parameter parameter_list) {
    // this function loads the hydro files
    hydro_type = paraRdr->getVal("hydro_type");
    load_viscous = paraRdr->getVal("load_viscous_info");
    if (hydro_type == 1) {
        HydroinfoH5* hydroinfo_ptr = new HydroinfoH5("JetData.h5", 500,
                                                     load_viscous);
        hydroinfo_ptr_in = hydroinfo_ptr;
    } else if (hydro_type == 2) {
        Hydroinfo_MUSIC* hydroinfo_ptr = new Hydroinfo_MUSIC();
        int hydro_mode = 8;
        int nskip_tau = paraRdr->getVal("hydro_nskip_tau");
        hydroinfo_ptr->readHydroData(hydro_mode, nskip_tau);
        hydroinfo_ptr_in = hydroinfo_ptr;
    } else if (hydro_type == 3) {
        Hydroinfo_MUSIC* hydroinfo_ptr = new Hydroinfo_MUSIC();
        int hydro_mode = 9;
        int nskip_tau = paraRdr->getVal("hydro_nskip_tau");
        hydroinfo_ptr->readHydroData(hydro_mode, nskip_tau);
        hydroinfo_ptr_in = hydroinfo_ptr;
    } else if (hydro_type == 4) {
        Hydroinfo_MUSIC* hydroinfo_ptr = new Hydroinfo_MUSIC();
        int hydro_mode = 10;
        int nskip_tau = 1;
        hydroinfo_ptr->readHydroData(hydro_mode, nskip_tau);
        hydroinfo_ptr_in = hydroinfo_ptr;
    } else {
        cout << "main: unrecognized hydro_type = " << hydro_type << endl;
        exit(1);
    }
    hydro_status = 1;
}


void HydroFile::evolve_hydro() {
    hydro_status = 3;
}

