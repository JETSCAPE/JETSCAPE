/////////////////////////////////////////////////////////////////////////
//                      hydrodynamics analysis
//
//              author: Chun Shen <chunshen@physics.mcgill.ca>
//              Copyright: Chun Shen 2014
//
//  This program load hydrodynamic evolution files and perform various
//  kinds of analysis
//  
//
//  To do in the future:
/////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>

//#include "./Hydroinfo_h5.h"
#include "./Hydroinfo_MUSIC.h"
#include "./Stopwatch.h"
#include "./FluidcellStatistic.h"
#include "./ParameterReader.h"
#include "./SurfaceFinder.h"

using namespace std;

int main(int argc, char *argv[]) {
    ParameterReader *paraRdr = new ParameterReader();
    paraRdr->readFromFile("parameters.dat");
    paraRdr->readFromArguments(argc, argv);
    paraRdr->echo();

    int load_viscous = paraRdr->getVal("load_viscous_info");
    int hydro_type = paraRdr->getVal("hydro_type");

    void* hydroinfo_ptr_in = NULL;

    Stopwatch sw;
    sw.tic();
    // hydro data file pointer
    if (hydro_type == 1) {
        //HydroinfoH5* hydroinfo_ptr = new HydroinfoH5("JetData.h5", 500,
        //                                             load_viscous);
        //hydroinfo_ptr_in = hydroinfo_ptr;
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

    FluidcellStatistic fluidcellanalysis(hydroinfo_ptr_in, paraRdr);
    double T_cut = paraRdr->getVal("T_cut");
    fluidcellanalysis.outputTempasTauvsX();
    fluidcellanalysis.outputKnudersonNumberasTauvsX();
    fluidcellanalysis.outputinverseReynoldsNumberasTauvsX();
    fluidcellanalysis.analysis_hydro_volume_for_photon(T_cut);
    fluidcellanalysis.output_temperature_vs_avg_utau();
    fluidcellanalysis.output_flowvelocity_vs_tau();

    // construct freeze-out hyper-surface
    // SurfaceFinder* surface_ptr = new SurfaceFinder(hydroinfo_ptr, paraRdr);
    // surface_ptr->Find_full_hypersurface();

    sw.toc();
    cout << "totally takes : " << sw.takeTime() << " seconds." << endl;

    return(0);
}


