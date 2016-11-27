// Copyright JETSCAPE Collaboration @ 2016
// This is a general basic class for hydrodynamics
// This is written by Longgang Pang and Chun Shen

#ifndef SRC_FLUID_DYNAMICS_H_
#define SRC_FLUID_DYNAMICS_H_

#include <vector>
#include <cstring>
#include <stdexcept>
#include <cmath>

#include "realtype.h"

enum HydroStatus {NOT_START, INITIALIZED, EVOLVING, FINISHED};

class JetSource {
    public:
        JetSource():j0(0.), j1(0.), j2(0.), j3(0.) {}
    private:
        real j0, j1, j2, j3;
};

//overload +-*/ for easier linear interpolation
class FluidCellInfo {
 public:
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

    FluidCellInfo() = default;    

    FluidCellInfo inline operator*=(real b);
};

/// adds \f$ c = a + b \f$
inline FluidCellInfo operator+(FluidCellInfo a, const FluidCellInfo & b) {
    a.energy_density += b.energy_density;
    a.entropy_density += b.entropy_density;
    a.temperature += b.temperature;
    a.pressure += b.pressure;
    a.qgp_fraction += b.qgp_fraction;
    a.mu_B += b.mu_B;
    a.mu_C += b.mu_C;
    a.mu_S += b.mu_S;
    a.vx += b.vx;
    a.vy += b.vy;
    a.vz += b.vz;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            a.pi[i][j] += b.pi[i][j];
        }
    }
    a.bulk_Pi += b.bulk_Pi;
    return a;
}

/// multiply the fluid cell with a scalar factor
FluidCellInfo inline FluidCellInfo::operator*=(real b){
    this->energy_density *= b;
    this->entropy_density *= b;
    this->temperature *= b;
    this->pressure *= b;
    this->qgp_fraction *= b;
    this->mu_B *= b;
    this->mu_C *= b;
    this->mu_S *= b;
    this->vx *= b;
    this->vy *= b;
    this->vz *= b;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            this->pi[i][j] *= b;
        }
    }
    this->bulk_Pi *= b;
    return *this;
}

/// multiply \f$ c = a * b \f$
inline FluidCellInfo operator*(real a, FluidCellInfo b){
    b *= a;
    return b;
}

/// multiply \f$ c = a * b \f$
inline FluidCellInfo operator*(FluidCellInfo a, real b){
    a *= b;
    return a;
}

/// division \f$ c = a / b \f$
inline FluidCellInfo operator/(FluidCellInfo a, real b){
    a *= 1.0/b;
    return a;
}

// print the fluid cell information for debuging
// this function has bugs
//std::ostream &operator<<(std::ostream &os, const FluidCellInfo &cell) {
//    os << "energy_density=" << cell.energy_density << std::endl; 
//    os << "entropy_density=" << cell.entropy_density << std::endl; 
//    os << "temperature=" << cell.temperature << std::endl; 
//    os << "pressure=" << cell.pressure << std::endl;
//    os << "qgp_fraction=" << cell.qgp_fraction << std::endl;
//    os << "mu_B=" << cell.mu_B << std::endl;
//    os << "mu_C=" << cell.mu_C << std::endl;
//    os << "mu_S=" << cell.mu_S << std::endl;
//    os << "vx=" << cell.vx << std::endl;
//    os << "vy=" << cell.vy << std::endl;
//    os << "vz=" << cell.vz << std::endl;
//    os << "pi[mu][nu]=" << std::endl;
//    for (int i = 0; i < 4; i++) {
//        for (int j = 0; j < 4; j++) {
//            os << cell.pi[i][j] << ' ';
//        }
//        os << std::endl;
//    }
//    os << "bulk_Pi=" << cell.bulk_Pi;
//    return os << std::endl;
//}

typedef struct {
    // data structure for outputing hyper-surface information
    real d3sigma_mu[4];     // surface vector
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
} SurfaceCellInfo;


class Parameter{
 public:
    // hydro parameters
    char* hydro_input_filename;
};


class EvolutionHistory{
 public:
    real tau_min, dtau;
    real x_min, dx;
    real y_min, dy;
    real eta_min, deta;
    int ntau, nx, ny, neta;
    // tau_eta_is_tz: default false;
    // true if hydro is in (t,x,y,z) coordinates
    bool tau_eta_is_tz;
    // the bulk information
    std::vector<FluidCellInfo> data;

    EvolutionHistory() {};
    ~EvolutionHistory() {data.clear();}

    inline real tau_max(){return tau_min + ntau * dtau;}
    inline real x_max(){return x_min + nx * dx;}
    inline real y_max(){return y_min + ny * dy;}
    inline real eta_max(){return eta_min + neta * deta;}

    class InvalidSpaceTimeRange : public std::invalid_argument {
        using std::invalid_argument::invalid_argument;
    };

    /** make sure the space time point (tau, x, y, eta) is inside
     * evolution history */
    void check_in_range(real tau, real x, real y, real eta) {
        if (tau < tau_min || tau > tau_max()) {
            throw InvalidSpaceTimeRange("tau=" + std::to_string(tau)
                    + " is not in range [" + std::to_string(tau_min) + "," 
                    + std::to_string(tau_max()) + "]");
        }
        if (x < x_min || x > x_max()) {
            throw InvalidSpaceTimeRange("x=" + std::to_string(x)
                    + " is not in range [" + std::to_string(x_min) + "," 
                    + std::to_string(x_max()) + "]");
        }
        if (y < y_min || y > y_max()) {
            throw InvalidSpaceTimeRange("y=" + std::to_string(y)
                    + " is not in range [" + std::to_string(y_min) + "," 
                    + std::to_string(y_max()) + "]");
        }
        if (eta < eta_min || eta > eta_max()) {
            throw InvalidSpaceTimeRange("eta=" + std::to_string(eta)
                    + " is not in range [" + std::to_string(eta_min) + "," 
                    + std::to_string(eta_max()) + "]");
        }
    }


    // get the lower bound of the fluid cell along tau
    inline int get_id_tau(real tau){
        return(static_cast<int>((tau - tau_min)/dtau));
    }
    // get the lower bound of the fluid cell along x
    inline int get_id_x(real x) { 
        return(static_cast<int>((x - x_min)/dx));
    }
    // get the lower bound of the fluid cell along y
    inline int get_id_y(real y) {
        return(static_cast<int>((y - y_min)/dy));
    }
    // get the lower bound of the fluid cell along y
    inline int get_id_eta(real eta) {
        return(static_cast<int>((eta - eta_min)/deta));
    }

    // get the coordinate of tau, x, y, eta on grid
    inline real tau_coord(int id_tau) { return tau_min + id_tau * dtau; }
    inline real x_coord(int id_x) { return x_min + id_x * dx; }
    inline real y_coord(int id_y) { return y_min + id_y * dy; }
    inline real eta_coord(int id_eta) { return eta_min + id_eta * deta; }

    // get the FluidCellInfo index in data
    inline int cell_index(int id_tau, int id_x, int id_y, int id_eta) {
        return  id_tau * nx * ny * neta + id_x * ny * neta
                        + id_y * neta + id_eta;
    }

    // get the FluidCellInfo at space point given time step
    FluidCellInfo get_at_time_step(int id_tau, real x, real y, real etas);

    // get the FluidCellInfo at given space time point
    FluidCellInfo get(real tau, real x, real y, real etas);
};


class FluidDynamics{
 protected:
    // record hydro start and end proper time [fm/c]
    real hydro_tau_0, hydro_tau_max;
    // record hydro freeze out temperature [GeV]
    real hydro_freeze_out_temperature;

    // record hydro running status
    HydroStatus hydro_status;  // 0: nothing happened
                       // 1: hydro has been initialized
                       // 2: hydro is evolving
                       // 3: all fluid cells have reached freeze-out, EvolutionHistory filled
                       // -1: An error occurred

 public:
    FluidDynamics() {};
    ~FluidDynamics() {};

    // How to store this data? In memory or hard disk?
    // 3D hydro may eat out the memory,
    // for large dataset, std::deque is better than std::vector.
    EvolutionHistory bulk_info;

    /*Keep this interface open in the beginning.*/
    // FreezeOutHyperSf hyper_sf;

    /* currently we have no standard for passing configurations */
    // pure virtual function; to be implemented by users
    // should make it easy to save evolution history to bulk_info
    virtual void initialize_hydro(Parameter parameter_list) {};

    virtual void evolve_hydro() {};

    virtual void evolve_hydro_one_step(JetSource jmu) {};

    // the following functions should be implemented in Jetscape
    int get_hydro_status() {return(hydro_status);}
    real get_hydro_start_time() {return(hydro_tau_0);}
    real get_hydro_end_time() {return(hydro_tau_max);}
    real get_hydro_freeze_out_temperature() {
        return(hydro_freeze_out_temperature);
    }

    /** retrive hydro information at a given space-tim point
     * throw InvalidSpaceTimeRange exception when
     * (t, x, y, z) is out of the EvolutionHistory range
     */
    virtual void get_hydro_info(real t, real x, real y, real z,
                                FluidCellInfo* fluid_cell_info_ptr){
        if (hydro_status != FINISHED || bulk_info.data.size() == 0) {
            throw std::runtime_error("Hydro evolution is not finished "
                                     "or EvolutionHistory is empty");
        }
        // judge whether to use 2D interpolation or 3D interpolation
        if (!bulk_info.tau_eta_is_tz) {
            real tau = std::sqrt(t * t - z * z);
            real eta = 0.5 * (std::log(t + z) - std::log(t - z));
            bulk_info.check_in_range(tau, x, y, eta);
            //return bulk_info.get(tau, x, y, eta);
        } else {
            bulk_info.check_in_range(t, x, y, z);
            //return bulk_info.get(t, x, y, z);
        }
    }

    // this function print out the information of the fluid cell to the screen
    void print_fluid_cell_information(FluidCellInfo* fluid_cell_info_ptr);

    // this function returns hypersurface for Cooper-Frye or recombination
    // the detailed implementation is left to the hydro developper
    virtual void get_hypersurface(real T_cut,
                                  SurfaceCellInfo* surface_list_ptr) {};

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

    real get_net_baryon_density(real time, real x, real y, real z);
    real get_net_charge_density(real time, real x, real y, real z);
};

#endif  // SRC_FLUID_DYNAMICS_H_
