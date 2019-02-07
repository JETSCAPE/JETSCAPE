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

#include "RealType.h"
#include "SurfaceFinder.h"
#include "cornelius.h"
#include "EvolutionHistory.h"

namespace Jetscape {

SurfaceFinder::SurfaceFinder(const Jetscape::real T_in,
                             const EvolutionHistory &bulk_data) :
    bulk_info(bulk_data) {

    T_cut = T_in;
}

/*
bool SurfaceFinder::check_intersect(
            Jetscape::real tau, Jetscape::real x, Jetscape::real y,
            Jetscape::real dt, Jetscape::real dx, Jetscape::real dy,
            Jetscape::real ***cube) {
    fluidCell *fluidCellptr = new fluidCell();
    bool intersect = true;

    auto tau_low = tau - dt/2.;
    auto tau_high = tau + dt/2.;
    auto x_left = x - dx/2.;
    auto x_right = x + dx/2.;
    auto y_left = y - dy/2.;
    auto y_right = y + dy/2.;

    if (hydro_type == 0) {
        hydroinfo_ptr->getHydroinfo(tau_low, x_left, y_left, fluidCellptr);
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_left, y_left, 0.0, tau_low,
                                            fluidCellptr);
    }
    cube[0][0][0] = fluidCellptr->temperature;
    if (hydro_type == 0) {
        hydroinfo_ptr->getHydroinfo(tau_low, x_left, y_right, fluidCellptr);
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_left, y_right, 0.0, tau_low,
                                            fluidCellptr);
    }
    cube[0][0][1] = fluidCellptr->temperature;
    if (hydro_type == 0) {
        hydroinfo_ptr->getHydroinfo(tau_low, x_right, y_left, fluidCellptr);
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_right, y_left, 0.0, tau_low,
                                            fluidCellptr);
    }
    cube[0][1][0] = fluidCellptr->temperature;
    if (hydro_type == 0) {
        hydroinfo_ptr->getHydroinfo(tau_low, x_right, y_right, fluidCellptr);
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_right, y_right, 0.0, tau_low,
                                            fluidCellptr);
    }
    cube[0][1][1] = fluidCellptr->temperature;
    if (hydro_type == 0) {
        hydroinfo_ptr->getHydroinfo(tau_high, x_left, y_left, fluidCellptr);
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_left, y_left, 0.0, tau_high,
                                            fluidCellptr);
    }
    cube[1][0][0] = fluidCellptr->temperature;
    if (hydro_type == 0) {
        hydroinfo_ptr->getHydroinfo(tau_high, x_left, y_right, fluidCellptr);
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_left, y_right, 0.0, tau_high,
                                            fluidCellptr);
    }
    cube[1][0][1] = fluidCellptr->temperature;
    if (hydro_type == 0) {
        hydroinfo_ptr->getHydroinfo(tau_high, x_right, y_left, fluidCellptr);
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_right, y_left, 0.0, tau_high,
                                            fluidCellptr);
    }
    cube[1][1][0] = fluidCellptr->temperature;
    if (hydro_type == 0) {
        hydroinfo_ptr->getHydroinfo(tau_high, x_right, y_right, fluidCellptr);
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_right, y_right, 0.0, tau_high,
                                            fluidCellptr);
    }
    cube[1][1][1] = fluidCellptr->temperature;

    if ((T_cut - cube[0][0][0])*(cube[1][1][1] - T_cut) < 0.0)
        if ((T_cut - cube[0][1][0])*(cube[1][0][1] - T_cut) < 0.0)
            if ((T_cut - cube[0][1][1])*(cube[1][0][0] - T_cut) < 0.0)
                if ((T_cut - cube[0][0][1])*(cube[1][1][0] - T_cut) < 0.0)
                    intersect = false;

    delete fluidCellptr;
    return(intersect);
}

int SurfaceFinder::Find_full_hypersurface() {
    ofstream output;
    output.open("hyper_surface_2+1d.dat");

    double grid_tau0, grid_tauf, grid_x0, grid_y0;
    if (hydro_type == 1) {
        grid_tau0 = hydroinfo_ptr->getHydrogridTau0();
        grid_tauf = hydroinfo_ptr->getHydrogridTaumax();
        grid_x0 = hydroinfo_ptr->getHydrogridX0();
        grid_y0 = hydroinfo_ptr->getHydrogridY0();
    } else {
        grid_tau0 = hydroinfo_MUSIC_ptr->get_hydro_tau0();
        grid_tauf = hydroinfo_MUSIC_ptr->get_hydro_tau_max();
        grid_x0 = (- hydroinfo_MUSIC_ptr->get_hydro_x_max()
                   + hydroinfo_MUSIC_ptr->get_hydro_dx());
        grid_y0 = grid_x0;
    }

    double grid_dt = paraRdr->getVal("grid_dt");
    double grid_dx = paraRdr->getVal("grid_dx");
    double grid_dy = paraRdr->getVal("grid_dy");

    int dim = 3;
    double *lattice_spacing = new double [dim];
    lattice_spacing[0] = grid_dt;
    lattice_spacing[1] = grid_dx;
    lattice_spacing[2] = grid_dy;

    Cornelius* cornelius_ptr = new Cornelius();
    cornelius_ptr->init(dim, T_cut, lattice_spacing);
  
    int ntime = static_cast<int>((grid_tauf - grid_tau0)/grid_dt);
    int nx = static_cast<int>(fabs(2.*grid_x0)/grid_dx);
    int ny = static_cast<int>(fabs(2.*grid_y0)/grid_dy);

    double ***cube = new double** [2];
    for (int i = 0; i < 2; i++) {
        cube[i] = new double* [2];
        for (int j = 0; j < 2; j++) {
            cube[i][j] = new double [2];
            for (int k = 0; k < 2; k++)
                cube[i][j][k] = 0.0;
        }
    }
    
    fluidCell *fluidCellptr = new fluidCell();
  
    for (int itime = 0; itime < ntime; itime++) {
        // loop over time evolution
        double tau_local = grid_tau0 + (itime + 0.5)*grid_dt;
        for (int i = 0; i < nx; i++) {
            // loops over the transverse plane
            double x_local = grid_x0 + (i + 0.5)*grid_dx;
            for (int j = 0; j < ny; j++) {
                double y_local = grid_y0 + (j + 0.5)*grid_dy;
                bool intersect = check_intersect(T_cut, tau_local, x_local,
                                                 y_local, grid_dt, grid_dx,
                                                 grid_dy, cube);
                if (intersect) {
                    cornelius_ptr->find_surface_3d(cube);
                    for (int isurf = 0; isurf < cornelius_ptr->get_Nelements();
                         isurf++) {
                        double tau_center = (
                                cornelius_ptr->get_centroid_elem(isurf, 0)
                                + tau_local - grid_dt/2.);
                        double x_center = (
                                cornelius_ptr->get_centroid_elem(isurf, 1)
                                + x_local - grid_dx/2.);
                        double y_center = (
                                cornelius_ptr->get_centroid_elem(isurf, 2)
                                + y_local - grid_dy/2.);

                        double da_tau =
                                cornelius_ptr->get_normal_elem(isurf, 0);
                        double da_x = cornelius_ptr->get_normal_elem(isurf, 1);
                        double da_y = cornelius_ptr->get_normal_elem(isurf, 2);
                       
                        if (hydro_type == 1) {
                            hydroinfo_ptr->getHydroinfo(
                                tau_center, x_center, y_center, fluidCellptr);
                        } else {
                            hydroinfo_MUSIC_ptr->getHydroValues(
                                x_center, y_center, 0.0, tau_center,
                                fluidCellptr);
                        }

                        output << scientific << setw(18) << setprecision(8) 
                               << tau_center << "   " << x_center << "   "
                               << y_center << "   " 
                               << da_tau << "   " << da_x << "   "
                               << da_y << "   " 
                               << fluidCellptr->temperature << "   "
                               << fluidCellptr->vx << "   " << fluidCellptr->vy 
                               << endl;

                    }
                }
            }
        }
    }
    output.close();
    
    delete fluidCellptr;
    delete cornelius_ptr;
    delete [] lattice_spacing;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++)
            delete [] cube[i][j];
        delete [] cube[i];
    }
    delete [] cube;
    return 0;
}
*/

}
