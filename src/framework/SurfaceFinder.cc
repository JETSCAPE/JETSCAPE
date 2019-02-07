/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/
// This is a general basic class for a hyper-surface finder

#include <cmath>
#include "RealType.h"
#include "SurfaceFinder.h"
#include "cornelius.h"
#include "EvolutionHistory.h"
#include "JetScapeLogger.h"

namespace Jetscape {

SurfaceFinder::SurfaceFinder(const Jetscape::real T_in,
                             const EvolutionHistory &bulk_data) :
    bulk_info(bulk_data) {

    T_cut = T_in;
    JSINFO << "Find a surface with temperature T = " << T_cut;
    boost_invariant = bulk_info.is_boost_invariant();
    if (boost_invariant) {
        JSINFO << "Hydro medium is boost invariant.";
    } else {
        JSINFO << "Hydro medium is not boost invariant.";
    }
    JSINFO << "Number of fluid cells = " << bulk_info.get_data_size();
}


SurfaceFinder::~SurfaceFinder() {
    surface_cell_list.clear();
}


void SurfaceFinder::Find_full_hypersurface() {
    if (boost_invariant) {
        JSINFO << "Finding a 2+1D hyper-surface at T = " << T_cut
               << " GeV ...";
        Find_full_hypersurface_3D();
    } else {
        JSINFO << "Finding a 3+1D hyper-surface at T = " << T_cut
               << " GeV ...";
        Find_full_hypersurface_4D();
    }
}


bool SurfaceFinder::check_intersect_3D(
            Jetscape::real tau, Jetscape::real x, Jetscape::real y,
            Jetscape::real dt, Jetscape::real dx, Jetscape::real dy,
            double ***cube) {
    bool intersect = true;

    auto tau_low = tau - dt/2.;
    auto tau_high = tau + dt/2.;
    auto x_left = x - dx/2.;
    auto x_right = x + dx/2.;
    auto y_left = y - dy/2.;
    auto y_right = y + dy/2.;

    auto fluid_cell = bulk_info.get(tau_low, x_left, y_left, 0.0);
    cube[0][0][0] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_low, x_left, y_right, 0.0);
    cube[0][0][1] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_low, x_right, y_left, 0.0);
    cube[0][1][0] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_low, x_right, y_right, 0.0);
    cube[0][1][1] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_high, x_left, y_right, 0.0);
    cube[1][0][0] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_high, x_left, y_right, 0.0);
    cube[1][0][1] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_high, x_right, y_left, 0.0);
    cube[1][1][0] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_high, x_right, y_right, 0.0);
    cube[1][1][1] = fluid_cell.temperature;

    if ((T_cut - cube[0][0][0])*(cube[1][1][1] - T_cut) < 0.0)
        if ((T_cut - cube[0][1][0])*(cube[1][0][1] - T_cut) < 0.0)
            if ((T_cut - cube[0][1][1])*(cube[1][0][0] - T_cut) < 0.0)
                if ((T_cut - cube[0][0][1])*(cube[1][1][0] - T_cut) < 0.0)
                    intersect = false;

    return(intersect);
}


void SurfaceFinder::Find_full_hypersurface_3D() {
    auto grid_tau0 = bulk_info.Tau0();
    auto grid_tauf = bulk_info.TauMax();
    auto grid_x0   = bulk_info.XMin();
    auto grid_y0   = bulk_info.YMin();;
    
    Jetscape::real grid_dt = 0.1;
    Jetscape::real grid_dx = 0.2;
    Jetscape::real grid_dy = 0.2;

    const int dim = 3;
    double lattice_spacing[dim];
    lattice_spacing[0] = grid_dt;
    lattice_spacing[1] = grid_dx;
    lattice_spacing[2] = grid_dy;

    std::unique_ptr<Cornelius> cornelius_ptr(new Cornelius());
    cornelius_ptr->init(dim, T_cut, lattice_spacing);
  
    const int ntime = static_cast<int>((grid_tauf - grid_tau0)/grid_dt);
    const int nx    = static_cast<int>(std::abs(2.*grid_x0)/grid_dx);
    const int ny    = static_cast<int>(std::abs(2.*grid_y0)/grid_dy);

    double ***cube = new double** [2];
    for (int i = 0; i < 2; i++) {
        cube[i] = new double* [2];
        for (int j = 0; j < 2; j++) {
            cube[i][j] = new double [2];
            for (int k = 0; k < 2; k++)
                cube[i][j][k] = 0.0;
        }
    }
    
    for (int itime = 0; itime < ntime; itime++) {
        // loop over time evolution
        auto tau_local = grid_tau0 + (itime + 0.5)*grid_dt;
        for (int i = 0; i < nx; i++) {
            // loops over the transverse plane
            auto x_local = grid_x0 + (i + 0.5)*grid_dx;
            for (int j = 0; j < ny; j++) {
                auto y_local = grid_y0 + (j + 0.5)*grid_dy;
                bool intersect = check_intersect_3D(tau_local, x_local,
                                                    y_local, grid_dt, grid_dx,
                                                    grid_dy, cube);
                if (intersect) {
                    cornelius_ptr->find_surface_3d(cube);
                    for (int isurf = 0; isurf < cornelius_ptr->get_Nelements();
                         isurf++) {
                        auto tau_center = (
                                cornelius_ptr->get_centroid_elem(isurf, 0)
                                + tau_local - grid_dt/2.);
                        auto x_center = (
                                cornelius_ptr->get_centroid_elem(isurf, 1)
                                + x_local - grid_dx/2.);
                        auto y_center = (
                                cornelius_ptr->get_centroid_elem(isurf, 2)
                                + y_local - grid_dy/2.);

                        auto da_tau =
                                cornelius_ptr->get_normal_elem(isurf, 0);
                        auto da_x = cornelius_ptr->get_normal_elem(isurf, 1);
                        auto da_y = cornelius_ptr->get_normal_elem(isurf, 2);
                       
                        auto fluid_cell = bulk_info.get(
                                        tau_center, x_center, y_center, 0.0);
                        auto surface_cell = PrepareASurfaceCell(
                                tau_center, x_center, y_center, 0.0,
                                da_tau, da_x, da_y, 0.0,
                                fluid_cell);
                        surface_cell_list.push_back(surface_cell);
                    }
                }
            }
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++)
            delete [] cube[i][j];
        delete [] cube[i];
    }
    delete [] cube;
}


bool SurfaceFinder::check_intersect_4D(
            Jetscape::real tau, Jetscape::real x, Jetscape::real y,
            Jetscape::real eta,
            Jetscape::real dt, Jetscape::real dx, Jetscape::real dy,
            Jetscape::real deta, double ****cube) {

    bool intersect = true;

    auto tau_low   = tau - dt/2.;
    auto tau_high  = tau + dt/2.;
    auto x_left    = x - dx/2.;
    auto x_right   = x + dx/2.;
    auto y_left    = y - dy/2.;
    auto y_right   = y + dy/2.;
    auto eta_left  = eta - deta/2.;
    auto eta_right = eta + deta/2.;

    auto fluid_cell = bulk_info.get(tau_low, x_left, y_left, eta_left);
    cube[0][0][0][0] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_low, x_left, y_left, eta_right);
    cube[0][0][0][1] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_low, x_left, y_right, eta_left);
    cube[0][0][1][0] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_low, x_left, y_right, eta_right);
    cube[0][0][1][1] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_low, x_right, y_left, eta_left);
    cube[0][1][0][0] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_low, x_right, y_left, eta_right);
    cube[0][1][0][1] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_low, x_right, y_right, eta_left);
    cube[0][1][1][0] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_low, x_right, y_right, eta_right);
    cube[0][1][1][1] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_high, x_left, y_left, eta_left);
    cube[1][0][0][0] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_high, x_left, y_left, eta_right);
    cube[1][0][0][1] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_high, x_left, y_right, eta_left);
    cube[1][0][1][0] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_high, x_left, y_right, eta_right);
    cube[1][0][1][1] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_high, x_right, y_left, eta_left);
    cube[1][1][0][0] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_high, x_right, y_left, eta_right);
    cube[1][1][0][1] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_high, x_right, y_right, eta_left);
    cube[1][1][1][0] = fluid_cell.temperature;
    fluid_cell = bulk_info.get(tau_high, x_right, y_right, eta_right);
    cube[1][1][1][1] = fluid_cell.temperature;

    if ((T_cut - cube[0][0][0][0])*(cube[1][1][1][1] - T_cut) < 0.0)
    if ((T_cut - cube[0][0][1][1])*(cube[1][1][0][0] - T_cut) < 0.0)
    if ((T_cut - cube[0][1][0][1])*(cube[1][0][1][0] - T_cut) < 0.0)
    if ((T_cut - cube[0][1][1][0])*(cube[1][0][0][1] - T_cut) < 0.0)
    if ((T_cut - cube[0][0][0][1])*(cube[1][1][1][0] - T_cut) < 0.0)
    if ((T_cut - cube[0][0][1][0])*(cube[1][1][0][1] - T_cut) < 0.0)
    if ((T_cut - cube[0][1][0][0])*(cube[1][0][1][1] - T_cut) < 0.0)
    if ((T_cut - cube[0][1][1][1])*(cube[1][0][0][0] - T_cut) < 0.0)
        intersect = false;

    return(intersect);
}

void SurfaceFinder::Find_full_hypersurface_4D() {
    auto grid_tau0 = bulk_info.Tau0();
    auto grid_tauf = bulk_info.TauMax();
    auto grid_x0   = bulk_info.XMin();
    auto grid_y0   = bulk_info.YMin();;
    auto grid_eta0 = bulk_info.EtaMin();;
    
    Jetscape::real grid_dt   = 0.1;
    Jetscape::real grid_dx   = 0.2;
    Jetscape::real grid_dy   = 0.2;
    Jetscape::real grid_deta = 0.2;

    const int dim = 4;
    double lattice_spacing[dim];
    lattice_spacing[0] = grid_dt;
    lattice_spacing[1] = grid_dx;
    lattice_spacing[2] = grid_dy;
    lattice_spacing[3] = grid_deta;

    std::unique_ptr<Cornelius> cornelius_ptr(new Cornelius());
    cornelius_ptr->init(dim, T_cut, lattice_spacing);
  
    const int ntime = static_cast<int>((grid_tauf - grid_tau0)/grid_dt);
    const int nx    = static_cast<int>(std::abs(2.*grid_x0)/grid_dx);
    const int ny    = static_cast<int>(std::abs(2.*grid_y0)/grid_dy);
    const int neta  = static_cast<int>(std::abs(2.*grid_eta0)/grid_deta);

    double ****cube = new double*** [2];
    for (int i = 0; i < 2; i++) {
        cube[i] = new double** [2];
        for (int j = 0; j < 2; j++) {
            cube[i][j] = new double* [2];
            for (int k = 0; k < 2; k++) {
                cube[i][j][k] = new double [2];
                for (int l = 0; l < 2; l++) {
                    cube[i][j][k][l] = 0.0;
                }
            }
        }
    }
    
    for (int itime = 0; itime < ntime; itime++) {
        // loop over time evolution
        auto tau_local = grid_tau0 + (itime + 0.5)*grid_dt;
        for (int l = 0; l < neta; l++) {
            auto eta_local = grid_eta0 + (l + 0.5)*grid_deta;
            // loops over the transverse plane
            for (int i = 0; i < nx; i++) {
                // loops over the transverse plane
                auto x_local = grid_x0 + (i + 0.5)*grid_dx;
                for (int j = 0; j < ny; j++) {
                    auto y_local = grid_y0 + (j + 0.5)*grid_dy;
                    bool intersect = check_intersect_4D(
                            tau_local, x_local, y_local, eta_local,
                            grid_dt, grid_dx, grid_dy, grid_deta, cube);
                    if (intersect) {
                        cornelius_ptr->find_surface_4d(cube);
                        for (int isurf = 0;
                             isurf < cornelius_ptr->get_Nelements(); isurf++) {
                            auto tau_center = (
                                    cornelius_ptr->get_centroid_elem(isurf, 0)
                                    + tau_local - grid_dt/2.);
                            auto x_center = (
                                    cornelius_ptr->get_centroid_elem(isurf, 1)
                                    + x_local - grid_dx/2.);
                            auto y_center = (
                                    cornelius_ptr->get_centroid_elem(isurf, 2)
                                    + y_local - grid_dy/2.);
                            auto eta_center = (
                                    cornelius_ptr->get_centroid_elem(isurf, 3)
                                    + eta_local - grid_deta/2.);

                            auto da_tau = (
                                    cornelius_ptr->get_normal_elem(isurf, 0));
                            auto da_x   = (
                                    cornelius_ptr->get_normal_elem(isurf, 1));
                            auto da_y   = (
                                    cornelius_ptr->get_normal_elem(isurf, 2));
                            auto da_eta = (
                                    cornelius_ptr->get_normal_elem(isurf, 3));
                           
                            auto fluid_cell = bulk_info.get(
                                tau_center, x_center, y_center, eta_center);
                            auto surface_cell = PrepareASurfaceCell(
                                    tau_center, x_center, y_center, eta_center,
                                    da_tau, da_x, da_y, da_eta,
                                    fluid_cell);
                            surface_cell_list.push_back(surface_cell);
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                delete [] cube[i][j][k];
            }
            delete [] cube[i][j];
        }
        delete [] cube[i];
    }
    delete [] cube;
}

SurfaceCellInfo SurfaceFinder::PrepareASurfaceCell(
    Jetscape::real tau, Jetscape::real x, Jetscape::real y, Jetscape::real eta,
    Jetscape::real da0, Jetscape::real da1, Jetscape::real da2,
    Jetscape::real da3, const FluidCellInfo fluid_cell) {

    SurfaceCellInfo temp_cell;
    temp_cell.tau           = tau;
    temp_cell.x             = x;
    temp_cell.y             = y;
    temp_cell.eta           = eta;
    temp_cell.d3sigma_mu[0] = da0;
    temp_cell.d3sigma_mu[1] = da1;
    temp_cell.d3sigma_mu[2] = da2;
    temp_cell.d3sigma_mu[3] = da3;

    temp_cell.energy_density  = fluid_cell.energy_density;
    temp_cell.entropy_density = fluid_cell.entropy_density;
    temp_cell.temperature     = fluid_cell.temperature;
    temp_cell.pressure        = fluid_cell.pressure;
    temp_cell.qgp_fraction    = fluid_cell.qgp_fraction;
    temp_cell.mu_B            = fluid_cell.mu_B;
    temp_cell.mu_C            = fluid_cell.mu_C;
    temp_cell.mu_S            = fluid_cell.mu_S;

    temp_cell.vx = fluid_cell.vx;
    temp_cell.vy = fluid_cell.vy;
    temp_cell.vz = fluid_cell.vz;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            temp_cell.pi[i][j] = fluid_cell.pi[i][j];
        }
    }
    temp_cell.bulk_Pi = fluid_cell.bulk_Pi;

    return(temp_cell);
}


}
