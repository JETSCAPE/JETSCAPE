// Copyright Chun Shen @ 2015
#ifndef SRC_SurfaceFinder_H_
#define SRC_SurfaceFinder_H_

#ifdef USE_HDF5
#include "./Hydroinfo_h5.h"
#endif

#include "./Hydroinfo_MUSIC.h"
#include "./ParameterReader.h"

using namespace std;

class SurfaceFinder {
 private:
    int hydro_type;
#ifdef USE_HDF5
    HydroinfoH5 *hydroinfo_ptr;
#endif
    Hydroinfo_MUSIC *hydroinfo_MUSIC_ptr;
    ParameterReader *paraRdr;
    double T_cut;

 public:
    SurfaceFinder(void* hydroinfo_ptr_in, ParameterReader* paraRdr_in);
    SurfaceFinder(void* hydroinfo_ptr_in, ParameterReader* paraRdr_in,
                  double T_cut_in);
    ~SurfaceFinder();

    bool check_intersect(double T_cut, double tau, double x, double y,
                         double dt, double dx, double dy, double ***cube);
    int Find_full_hypersurface();
};

#endif  // SRC_SurfaceFinder_H_
