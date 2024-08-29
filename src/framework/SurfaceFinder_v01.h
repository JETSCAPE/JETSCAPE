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

#ifndef SurfaceFinder_v01_H_
#define SurfaceFinder_v01_H_

#include <vector>

#include "RealType.h"
#include "SurfaceCellInfo.h"
#include "FluidEvolutionHistory.h"

namespace Jetscape {

class SurfaceFinder_v01 {
private:
  Jetscape::real T_cut;
  const EvolutionHistory &bulk_info;
  bool boost_invariant;

  std::vector<SurfaceCellInfo> surface_cell_list;

public:
  SurfaceFinder_v01(const Jetscape::real T_in, const EvolutionHistory &bulk_data);
  ~SurfaceFinder_v01();

  void Find_full_hypersurface();

  int get_number_of_surface_cells() const { return (surface_cell_list.size()); }
  SurfaceCellInfo get_surface_cell_with_idx(int idx) const {
    return (surface_cell_list[idx]);
  }
  std::vector<SurfaceCellInfo> get_surface_cells_vector() const {
    return (surface_cell_list);
  }

#pragma region check_intersect_3D
/**
 * @brief Checks if there is an intersection in a 3D hypersurface.
 * 
 * @param tau Local tau value.
 * @param x Local x value.
 * @param y Local y value.
 * @param dt Time step.
 * @param dx X step.
 * @param dy Y step.
 * @param cube Temperature values at the corners of the cube.
 * @return True if there is an intersection, false otherwise.
 */
bool SurfaceFinder_v01::check_intersect_3D(Jetscape::real tau, Jetscape::real x,
                                       Jetscape::real y, Jetscape::real dt,
                                       Jetscape::real dx, Jetscape::real dy,
                                       std::array<std::array<std::array<double, 2>, 2>, 2>& cube);
/**
 * @brief Fills the 4D array cube with temperature values from the fluid cells.
 *
 * @param tau Central value of tau.
 * @param x Central value of x.
 * @param y Central value of y.
 * @param dt Time step size.
 * @param dx X step size.
 * @param dy Y step size.
 * @param cube 3D array to store temperature values of the grid cell.
 */
void SurfaceFinder_v01::fill_cube_with_temperatures(
  Jetscape::real tau, Jetscape::real x, Jetscape::real y, 
  Jetscape::real dt, Jetscape::real dx, Jetscape::real dy,
  std::array<std::array<std::array<double, 2>, 2>, 2>& cube);

  
/**
 * @brief Checks if the temperature values in the cube intersect the cutoff temperature.
 * 
 * @param cube Temperature values at the corners of the cube.
 * @return True if the temperature values intersect the cutoff temperature, false otherwise.
 */
bool SurfaceFinder_v01::temperature_intersects_cutoff(const std::array<std::array<std::array<double, 2>, 2>, 2>& cube);
#pragma endregion check_intersect_3D

#pragma region Find_fill_hypersurface_3D
/**
 * @brief Find the full hypersurface in 3D.
 */
void SurfaceFinder_v01::Find_full_hypersurface_3D();

/**
 * @brief Process a single time step in the hypersurface finding.
 *
 * @param itime Current time step index.
 * @param grid_tau0 Initial tau value.
 * @param grid_dt Tau step size.
 * @param nx Number of x steps.
 * @param ny Number of y steps.
 * @param grid_x0 Initial x value.
 * @param grid_dx X step size.
 * @param grid_y0 Initial y value.
 * @param grid_dy Y step size.
 * @param cube 3D array to store temperature values.
 * @param cornelius_ptr Pointer to Cornelius object.
 */
void SurfaceFinder_v01::process_time_step(int itime, Jetscape::real grid_tau0, Jetscape::real grid_dt, int nx, int ny, 
                                      Jetscape::real grid_x0, Jetscape::real grid_dx, Jetscape::real grid_y0, Jetscape::real grid_dy, 
                                      std::array<std::array<std::array<double, 2>, 2>, 2>& cube, 
                                      const std::unique_ptr<Cornelius>& cornelius_ptr);

/**
 * @brief Process the x-plane in the hypersurface finding.
 *
 * @param i Current x index.
 * @param tau_local Local tau value.
 * @param grid_x0 Initial x value.
 * @param grid_dx X step size.
 * @param ny Number of y steps.
 * @param grid_y0 Initial y value.
 * @param grid_dy Y step size.
 * @param cube 3D array to store temperature values.
 * @param cornelius_ptr Pointer to Cornelius object.
 */
void SurfaceFinder_v01::process_x_plane(int i, Jetscape::real tau_local, Jetscape::real grid_x0, Jetscape::real grid_dx, int ny, 
                                    Jetscape::real grid_y0, Jetscape::real grid_dy, 
                                    std::array<std::array<std::array<double, 2>, 2>, 2>& cube, 
                                    const std::unique_ptr<Cornelius>& cornelius_ptr);

/**
 * @brief Process the y-plane in the hypersurface finding.
 *
 * @param j Current y index.
 * @param tau_local Local tau value.
 * @param x_local Local x value.
 * @param grid_y0 Initial y value.
 * @param grid_dy Y step size.
 * @param cube 3D array to store temperature values.
 * @param cornelius_ptr Pointer to Cornelius object.
 */
void SurfaceFinder_v01::process_y_plane(int j, Jetscape::real tau_local, Jetscape::real x_local, Jetscape::real grid_y0, Jetscape::real grid_dy, 
                                    std::array<std::array<std::array<double, 2>, 2>, 2>& cube, 
                                    const std::unique_ptr<Cornelius>& cornelius_ptr);

  /**
 * @brief Process surface elements found by Cornelius.
 *
 * @param tau_local Local tau value.
 * @param x_local Local x value.
 * @param y_local Local y value.
 * @param grid_dt Tau step size.
 * @param grid_dx X step size.
 * @param grid_dy Y step size.
 * @param cube 3D array to store temperature values.
 * @param cornelius_ptr Pointer to Cornelius object.
 */
void SurfaceFinder_v01::process_surface_elements(Jetscape::real tau_local, Jetscape::real x_local, Jetscape::real y_local, 
                                             Jetscape::real grid_dt, Jetscape::real grid_dx, Jetscape::real grid_dy, 
                                             std::array<std::array<std::array<double, 2>, 2>, 2>& cube, 
                                             const std::unique_ptr<Cornelius>& cornelius_ptr);
                                            
#pragma endregion Find_fill_hypersurface_3D

#pragma region check_intersect_4D
/**
 * @brief Checks if the temperature in a 4D grid cell intersects a given temperature cutoff.
 *
 * @param tau Central value of tau.
 * @param x Central value of x.
 * @param y Central value of y.
 * @param eta Central value of eta.
 * @param dt Time step size.
 * @param dx X step size.
 * @param dy Y step size.
 * @param deta Eta step size.
 * @param cube 4D array to store temperature values of the grid cell.
 * @return True if the temperature intersects the cutoff, false otherwise.
 */
bool SurfaceFinder_v01::check_intersect_4D(
    Jetscape::real tau, Jetscape::real x, Jetscape::real y, Jetscape::real eta,
    Jetscape::real dt, Jetscape::real dx, Jetscape::real dy, Jetscape::real deta,
    std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> &cube);

/**
 * @brief Fills the 4D array cube with temperature values from the fluid cells.
 *
 * @param tau Central value of tau.
 * @param x Central value of x.
 * @param y Central value of y.
 * @param eta Central value of eta.
 * @param dt Time step size.
 * @param dx X step size.
 * @param dy Y step size.
 * @param deta Eta step size.
 * @param cube 4D array to store temperature values of the grid cell.
 */
void SurfaceFinder_v01::fill_cube_with_temperatures(
    Jetscape::real tau, Jetscape::real x, Jetscape::real y, Jetscape::real eta,
    Jetscape::real dt, Jetscape::real dx, Jetscape::real dy, Jetscape::real deta,
    std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> &cube);


/**
 * @brief Checks if the temperature in the 4D array cube intersects the cutoff temperature.
 *
 * @param cube 4D array containing temperature values of the grid cell.
 * @return True if the temperature intersects the cutoff, false otherwise.
 */
bool SurfaceFinder_v01::temperature_intersects_cutoff(
    const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> &cube);
    
#pragma endregion check_intersect_4D
  
#pragma region Find_full_hypersurface_4D
/**
 * @brief Finds the full hypersurface in 4D space.
 *
 * This method calculates the hypersurface by iterating over time, transverse plane, and rapidity
 * using a grid. It initializes the Cornelius object and processes each grid point to check if it
 * intersects with the hypersurface. If so, it computes the surface elements and stores them in
 * the `surface_cell_list`.
 */
void SurfaceFinder_v01::Find_full_hypersurface_4D();

/**
 * @brief Creates a 4D cube with initial values set to 0.0.
 *
 * @return A 4D vector representing the cube with dimensions [2][2][2][2] initialized to 0.0.
 */
std::vector<std::vector<std::vector<std::vector<double>>>> SurfaceFinder_v01::create_cube();

/**
 * @brief Processes a time slice of the hypersurface.
 *
 * Iterates over rapidity and transverse plane for a given time slice to check for intersections
 * and update surface cells.
 *
 * @param itime The current time index in the grid.
 * @param grid_tau0 The initial tau value of the grid.
 * @param grid_dt The time step size.
 * @param grid_eta0 The initial eta value of the grid.
 * @param grid_deta The eta step size.
 * @param neta The number of eta grid points.
 * @param grid_x0 The initial x value of the grid.
 * @param grid_dx The x step size.
 * @param nx The number of x grid points.
 * @param grid_y0 The initial y value of the grid.
 * @param grid_dy The y step size.
 * @param ny The number of y grid points.
 * @param cube The 4D vector representing the grid cube.
 * @param cornelius_ptr A unique pointer to the Cornelius object used for surface finding.
 */
void SurfaceFinder_v01::process_time_slice(
    int itime, Jetscape::real grid_tau0, Jetscape::real grid_dt,
    Jetscape::real grid_eta0, Jetscape::real grid_deta, int neta,
    Jetscape::real grid_x0, Jetscape::real grid_dx, int nx,
    Jetscape::real grid_y0, Jetscape::real grid_dy, int ny,
    std::vector<std::vector<std::vector<std::vector<double>>>>& cube,
    const std::unique_ptr<Cornelius>& cornelius_ptr);

/**
 * @brief Processes an eta slice of the hypersurface.
 *
 * Iterates over the transverse plane for a given eta slice to check for intersections
 * and update surface cells.
 *
 * @param l The current eta index in the grid.
 * @param tau_local The local tau value.
 * @param grid_eta0 The initial eta value of the grid.
 * @param grid_deta The eta step size.
 * @param grid_x0 The initial x value of the grid.
 * @param grid_dx The x step size.
 * @param nx The number of x grid points.
 * @param grid_y0 The initial y value of the grid.
 * @param grid_dy The y step size.
 * @param ny The number of y grid points.
 * @param cube The 4D vector representing the grid cube.
 * @param cornelius_ptr A unique pointer to the Cornelius object used for surface finding.
 */
void SurfaceFinder_v01::process_eta_slice(
    int l, Jetscape::real tau_local, Jetscape::real grid_eta0, Jetscape::real grid_deta,
    Jetscape::real grid_x0, Jetscape::real grid_dx, int nx,
    Jetscape::real grid_y0, Jetscape::real grid_dy, int ny,
    std::vector<std::vector<std::vector<std::vector<double>>>>& cube,
    const std::unique_ptr<Cornelius>& cornelius_ptr);

/**
 * @brief Processes an x slice of the hypersurface.
 *
 * Iterates over the y grid points for a given x slice to check for intersections
 * and update surface cells.
 *
 * @param i The current x index in the grid.
 * @param tau_local The local tau value.
 * @param eta_local The local eta value.
 * @param grid_x0 The initial x value of the grid.
 * @param grid_dx The x step size.
 * @param grid_y0 The initial y value of the grid.
 * @param grid_dy The y step size.
 * @param ny The number of y grid points.
 * @param cube The 4D vector representing the grid cube.
 * @param cornelius_ptr A unique pointer to the Cornelius object used for surface finding.
 */
void SurfaceFinder_v01::process_x_slice(
    int i, Jetscape::real tau_local, Jetscape::real eta_local, Jetscape::real grid_x0, Jetscape::real grid_dx,
    Jetscape::real grid_y0, Jetscape::real grid_dy, int ny,
    std::vector<std::vector<std::vector<std::vector<double>>>>& cube,
    const std::unique_ptr<Cornelius>& cornelius_ptr);

/**
 * @brief Processes a single grid point to check for intersection and compute surface cells.
 *
 * Checks if the grid point intersects with the hypersurface, and if so, computes the surface
 * elements and adds them to the `surface_cell_list`.
 *
 * @param tau_local The local tau value.
 * @param x_local The local x value.
 * @param y_local The local y value.
 * @param eta_local The local eta value.
 * @param grid_dt The time step size.
 * @param grid_dx The x step size.
 * @param grid_dy The y step size.
 * @param grid_deta The eta step size.
 * @param cube The 4D vector representing the grid cube.
 * @param cornelius_ptr A unique pointer to the Cornelius object used for surface finding.
 */
void SurfaceFinder_v01::process_grid_point(
    Jetscape::real tau_local, Jetscape::real x_local, Jetscape::real y_local, Jetscape::real eta_local,
    Jetscape::real grid_dt, Jetscape::real grid_dx, Jetscape::real grid_dy, Jetscape::real grid_deta,
    std::vector<std::vector<std::vector<std::vector<double>>>>& cube,
    const std::unique_ptr<Cornelius>& cornelius_ptr);

/**
 * @brief Computes the centroid positions of the surface elements.
 *
 * @param cornelius_ptr A unique pointer to the Cornelius object used for surface finding.
 * @param isurf The index of the surface element.
 * @param tau_local The local tau value.
 * @param x_local The local x value.
 * @param y_local The local y value.
 * @param eta_local The local eta value.
 * @param grid_dt The time step size.
 * @param grid_dx The x step size.
 * @param grid_dy The y step size.
 * @param grid_deta The eta step size.
 * @return A tuple containing the centroid positions (tau_center, x_center, y_center, eta_center).
 */
std::tuple<Jetscape::real, Jetscape::real, Jetscape::real, Jetscape::real> SurfaceFinder_v01::compute_centroids(
    const std::unique_ptr<Cornelius>& cornelius_ptr, int isurf, Jetscape::real tau_local, Jetscape::real x_local, Jetscape::real y_local, Jetscape::real eta_local,
    Jetscape::real grid_dt, Jetscape::real grid_dx, Jetscape::real grid_dy, Jetscape::real grid_deta);

/**
 * @brief Computes the normal vectors of the surface elements.
 *
 * @param cornelius_ptr A unique pointer to the Cornelius object used for surface finding.
 * @param isurf The index of the surface element.
 * @return A tuple containing the normal vectors (da_tau, da_x, da_y, da_eta).
 */
std::tuple<Jetscape::real, Jetscape::real, Jetscape::real, Jetscape::real> SurfaceFinder_v01::compute_normals(
    const std::unique_ptr<Cornelius>& cornelius_ptr, int isurf);
#pragma endregion Find_full_hypersurface_4D

/**
 * @brief Prepares a surface cell from given parameters and fluid cell information.
 *
 * This function initializes a `SurfaceCellInfo` structure using the provided
 * parameters and fluid cell data. It computes the four-velocity components 
 * and copies the stress-energy tensor components.
 *
 * @param tau The proper time.
 * @param x The x-coordinate.
 * @param y The y-coordinate.
 * @param eta The space-time rapidity.
 * @param da0 The 0th component of the surface normal vector.
 * @param da1 The 1st component of the surface normal vector.
 * @param da2 The 2nd component of the surface normal vector.
 * @param da3 The 3rd component of the surface normal vector.
 * @param fluid_cell The fluid cell information.
 * @return A populated `SurfaceCellInfo` structure.
 */
  SurfaceCellInfo PrepareASurfaceCell(Jetscape::real tau, Jetscape::real x,
                                      Jetscape::real y, Jetscape::real eta,
                                      Jetscape::real da0, Jetscape::real da1,
                                      Jetscape::real da2, Jetscape::real da3,
                                      const FluidCellInfo fluid_cell);
};
/**
 * @brief Populate the pi tensor for the surface cell.
 * 
 * @param temp_cell Reference to the surface cell.
 * @param fluid_cell Reference to the fluid cell.
 */
void SurfaceFinder_v01::populate_pi_tensor(SurfaceCellInfo& temp_cell, const FluidCellInfo& fluid_cell) ;
} // namespace Jetscape

#endif // SurfaceFinder_v01_H_
