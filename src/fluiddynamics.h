// Copyright JETSCAPE @ 2016

#ifndef FLUIDDYNAMICS_H_
#define FLUIDDYNAMICS_H_

#include <tuple>
#include <vector>

typedef float real

typedef std::tuple<real, real, real> real3;
typedef std::tuple<real, real, real, real> real4;

typedef struct {
    real energy_density;
    real temperature;
    real entropy_density;
    real qgp_fraction;
    real vx, vy, vz;
    // do we need pi^{mu nu}, Pi, net_baryon, net_charge?
    // for thermal photon or other kinds of studies
} BulkElement;

typedef struct {
    // data structure for outputing fluid cell information
    real energy_density;    // local energy density [GeV/fm^3]
    real entropy_density;   // local entropy density [1/fm^3]
    real temperature;       // local temperature [GeV]
    real pressure;          // thermal pressure [GeV/fm^3]
    real qgp_fraction;
    real mu_B;              // net baryon chemical potential [GeV]
    real mu_C;              // net charge chemical potential [GeV]
    real mu_S;              // net strangeness chemical potential [GeV]
    real vx, vy, vz;        // flow velocity
    real pi[4][4];          // shear stress tensor [GeV/fm^3]
    real bulk_Pi;           // bulk viscous pressure [GeV/fm^3]
} FluidCellInfo;

class EvolutionHistory{
 public:
    real tmin, nt, dt;
    real xmin, nx, dx;
    real ymin, ny, dy;
    real zmin, nz, dz;
    // default: set using_tz_for_tau_eta=true
    bool  using_tz_for_tau_eta;
    // the bulk information
    std::vector<BulkElement> data;

    EvolutionHistory();
};


class FluidDynamics{
 public:
    // How to store this data? In memory or hard disk?
    // 3D hydro may eat out the memory,
    // for large dataset, std::deque is better than std::vector.
    EvolutionHistory bulk_info;

    /*Keep this interface open in the beginning.*/
    // FreezeOutHyperSf hyper_sf;

    /* currently we have no standard for passing configurations */
    // pure virtual function; to be implemented by users
    // should make it easy to save evolution history to bulk_info
    virtual void evolution(const EnergyMomentumTensor & tmn,
                           real xmax, real ymax, real hmax,
                           real tau0, real dtau, real dx,
                           real dy, real dh, real etaos,
                           real dtau_out, real dx_out, real dy_out,
                           real dh_out) const = 0;

    // the following functions should be implemented in Jetscape

    // main function to retrive hydro information
    void get_hydro_info(real time, real x, real y, real z,
                        FluidCellInfo* fluid_cell_info_ptr);
    // all the following functions will call function get_hydro_info()
    // to get thermaldynamic and dynamical information at a space-time point
    // (time, x, y, z)

    real get_energy_density(real time, real x, real y, real z);
    real get_entropy_density(real time, real x, real y, real z);
    real get_temperature(real time, real x, real y, real z);
    real get_qgp_fraction(real time, real x, real y, real z);
    // real3 return std::make_tuple(vx, vy, vz)
    real3 get_3fluid_velocity(real time, real x, real y, real z);
    // real4 return std::make_tuple(ut, ux, uy, uz)
    real4 get_4fluid_velocity(real time, real x, real y, real z);
    // real get_net_baryon_density(real time, real x, real y, real z);
    // real get_net_charge_density(real time, real x, real y, real z);
};

#endif  // FLUIDDYNAMICS_H_
