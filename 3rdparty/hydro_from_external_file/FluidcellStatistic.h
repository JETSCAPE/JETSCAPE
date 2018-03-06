// Copyright Chun Shen @ 2015
#ifndef SRC_FLUIDCELLSTATISTIC_H_
#define SRC_FLUIDCELLSTATISTIC_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>

#ifdef USE_HDF5
#include "./Hydroinfo_h5.h"
#endif

#include "./Hydroinfo_MUSIC.h"
#include "./ParameterReader.h"

// using namespace std;
using std::string;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::abs;

class FluidcellStatistic {
 private:
    int hydro_type;
#ifdef USE_HDF5
    HydroinfoH5 *hydroinfo_ptr;
#endif
    Hydroinfo_MUSIC *hydroinfo_MUSIC_ptr;
    ParameterReader *paraRdr;
    double T_dec;
    double hbarC;
    double grid_dt, grid_dx, grid_dy;
    double grid_tau0, grid_tauf;
    double grid_x0, grid_y0;

 public:
    FluidcellStatistic(void* hydroinfo_ptr_in, ParameterReader* paraRdr_in);
    ~FluidcellStatistic();
    void checkFreezeoutSurface(double Tdec);
    double compute_local_expansion_rate(double tau_local, double x_local,
                                        double y_local);
    void output_momentum_anisotropy_vs_tau();
    void output_temperature_vs_tau();
    void output_flowvelocity_vs_tau();
    void output_temperature_vs_avg_utau();
    void outputTempasTauvsX();
    void outputKnudersonNumberasTauvsX();
    void outputinverseReynoldsNumberasTauvsX();
    void outputBulkinverseReynoldsNumberasTauvsX();
    void analysis_hydro_volume_for_photon(double T_cut);
    double calculate_spacetime_4volume(double T_cut);
    double calculate_average_tau(double T_cut);
    double calculate_average_temperature4(double T_cut);
    double calculate_average_integrated_photonRate_parameterization(
                                                                double T_cut);
    double calculate_hypersurface_3volume(double T_cut);
};

#endif  // SRC_FLUIDCELLSTATISTIC_H_
