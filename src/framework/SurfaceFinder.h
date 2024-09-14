/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
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

#ifndef SURFACEFINDER_H_
#define SURFACEFINDER_H_

#include <fstream>
#include <omp.h>
#include <vector>
#include <memory>

#include "FluidEvolutionHistory.h"
#include "RealType.h"
#include "SurfaceCellInfo.h"
#include "cornelius.h"

namespace Jetscape {

class SurfaceFinder {
 private:
  Jetscape::real T_cut;
  const EvolutionHistory &bulk_info;
  bool boost_invariant;

  std::vector<SurfaceCellInfo> surface_cell_list;

 public:
  SurfaceFinder(const Jetscape::real T_in, const EvolutionHistory &bulk_data);
  ~SurfaceFinder();

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
 * @brief Checks if the cube intersects the cutoff temperature.
 * 
 * @param tau Central value of tau.
 * @param x Central value of x.
 * @param y Central value of y.
 * @param dt Time step size.
 * @param cube 3D array to store temperature values of the grid cell.
 * @return True if the cube intersects the cutoff temperature, false otherwise.
 */
bool check_intersect_3D(Jetscape::real tau, Jetscape::real x,
                        Jetscape::real y, Jetscape::real dt,
                        Jetscape::real dx, Jetscape::real dy,
                        // std::array<std::array<std::array<double, 2>, 2>, 2>& cube
                        double ***cube
                        );

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
void fill_cube_with_temperatures(
  Jetscape::real tau, Jetscape::real x, Jetscape::real y, 
  Jetscape::real dt, Jetscape::real dx, Jetscape::real dy,
  // std::array<std::array<std::array<double, 2>, 2>, 2>& cube
  double ***cube
  );

/**
 * @brief Checks if the temperature values in the cube intersect the cutoff temperature.
 * 
 * @param cube Temperature values at the corners of the cube.
 * @return True if the temperature values intersect the cutoff temperature, false otherwise.
 */
bool temperature_intersects_cutoff(
  //const 
// std::array<std::array<std::array<double, 2>, 2>, 2>& cube
double ***cube
);
#pragma endregion check_intersect_3D
#pragma region Find_fill_hypersurface_3D
  void Find_full_hypersurface_3D();
  /*
* @brief Reduce the local 2D vector to 1D class member vector.
*
* @param surface_cell_list_local Local 2D vector.
* @param surface_cell_list_sz Size of the local 2D vector.
*/
void reduce_surface_cell_list(std::vector<std::vector<SurfaceCellInfo>>& surface_cell_list_local, int surface_cell_list_sz); 
  void process_surface_elements(Jetscape::real tau_local, Jetscape::real x_local, Jetscape::real y_local, 
                                             Jetscape::real grid_dt, Jetscape::real grid_dx, Jetscape::real grid_dy, 
                                            //  std::array<std::array<std::array<double, 2>, 2>, 2>& cube,
                                              double ***cube,
                                             const int itime, const int nx ,const int ny ,int i ,int j, 
                                             const std::unique_ptr<Cornelius>& cornelius_ptr, std::vector<std::vector<SurfaceCellInfo>>& surface_cell_list_local);
#pragma endregion Find_fill_hypersurface_3D
  bool check_intersect_4D(Jetscape::real tau, Jetscape::real x,
                          Jetscape::real y, Jetscape::real eta,
                          Jetscape::real dt, Jetscape::real dx,
                          Jetscape::real dy, Jetscape::real deta,
                          double ****cube);
  void Find_full_hypersurface_4D();

  SurfaceCellInfo PrepareASurfaceCell(Jetscape::real tau, Jetscape::real x,
                                      Jetscape::real y, Jetscape::real eta,
                                      Jetscape::real da0, Jetscape::real da1,
                                      Jetscape::real da2, Jetscape::real da3,
                                      const FluidCellInfo fluid_cell);

  void WriteSurfaceToFile(const std::vector<SurfaceCellInfo> &surface_cell_list,
                          std::string filename);
};

}  // namespace Jetscape

#endif  // SURFACEFINDER_H_
