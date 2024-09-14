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

#include "SurfaceFinder.h"

#include <cmath>

#include "FluidEvolutionHistory.h"
#include "JetScapeLogger.h"
#include "RealType.h"
#include "cornelius.h"

namespace Jetscape {

SurfaceFinder::SurfaceFinder(const Jetscape::real T_in,
                             const EvolutionHistory &bulk_data)
    : bulk_info(bulk_data) {
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

SurfaceFinder::~SurfaceFinder() { surface_cell_list.clear(); }

void SurfaceFinder::Find_full_hypersurface() {
  char* surf_path = std::getenv("SURF_PATH");
  std::string surf_path_str;
  std::string filename;
  if (surf_path != NULL) {
    JSINFO << "SURF_PATH is set to " << surf_path;
    surf_path_str = std::string(surf_path);
    filename = surf_path_str + "/hypersurface.dat";
  } else {
    JSINFO << "SURF_PATH is not set.";
    filename = "hypersurface.dat";
  }
  if (boost_invariant) {
    JSINFO << "Finding a 2+1D hyper-surface at T = " << T_cut << " GeV ...";
    auto start = std::chrono::high_resolution_clock::now();
    Find_full_hypersurface_3D();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    JSINFO << "3D Time to find the hypersurface: " << elapsed_seconds.count()
           << " s";
    WriteSurfaceToFile(surface_cell_list, filename);
  } else {
    JSINFO << "Finding a 3+1D hyper-surface at T = " << T_cut << " GeV ...";
    auto start = std::chrono::high_resolution_clock::now();
    Find_full_hypersurface_4D();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    JSINFO << "4D Time to find the hypersurface: " << elapsed_seconds.count()
           << " s";
    WriteSurfaceToFile(surface_cell_list, filename);
  }
}
#pragma region check_intersect_3D
/**
 * @brief Checks if the temperature values in the cube intersect the cutoff temperature.
 * 
 * @param tau Central value of tau.
 * @param x Central value of x.
 * @param y Central value of y.
 * @param dt Time step size.
 * @param dx X step size.
 * @param dy Y step size.
 * @param cube 3D array to store temperature values of the grid cell.
 * @return True if the temperature values intersect the cutoff temperature, false otherwise.
 */
bool SurfaceFinder::check_intersect_3D(Jetscape::real tau, Jetscape::real x,
                                       Jetscape::real y, Jetscape::real dt,
                                       Jetscape::real dx, Jetscape::real dy,
                                      //  std::array<std::array<std::array<double, 2>, 2>, 2>& cube
                                       double ***cube
                                       ) {

  bool intersect = true;

  auto tau_low = tau - dt / 2.;
  auto tau_high = tau + dt / 2.;
  auto x_left = x - dx / 2.;
  auto x_right = x + dx / 2.;
  auto y_left = y - dy / 2.;
  auto y_right = y + dy / 2.;

  auto fluid_cell = bulk_info.get(tau_low, x_left, y_left, 0.0);
  cube[0][0][0] = fluid_cell.temperature;
  fluid_cell = bulk_info.get(tau_low, x_left, y_right, 0.0);
  cube[0][0][1] = fluid_cell.temperature;
  fluid_cell = bulk_info.get(tau_low, x_right, y_left, 0.0);
  cube[0][1][0] = fluid_cell.temperature;
  fluid_cell = bulk_info.get(tau_low, x_right, y_right, 0.0);
  cube[0][1][1] = fluid_cell.temperature;
  fluid_cell = bulk_info.get(tau_high, x_left, y_left, 0.0);
  cube[1][0][0] = fluid_cell.temperature;
  fluid_cell = bulk_info.get(tau_high, x_left, y_right, 0.0);
  cube[1][0][1] = fluid_cell.temperature;
  fluid_cell = bulk_info.get(tau_high, x_right, y_left, 0.0);
  cube[1][1][0] = fluid_cell.temperature;
  fluid_cell = bulk_info.get(tau_high, x_right, y_right, 0.0);
  cube[1][1][1] = fluid_cell.temperature;

  if ((T_cut - cube[0][0][0]) * (cube[1][1][1] - T_cut) < 0.0)
    if ((T_cut - cube[0][1][0]) * (cube[1][0][1] - T_cut) < 0.0)
      if ((T_cut - cube[0][1][1]) * (cube[1][0][0] - T_cut) < 0.0)
        if ((T_cut - cube[0][0][1]) * (cube[1][1][0] - T_cut) < 0.0)
          intersect = false;

  return (intersect);
  // fill_cube_with_temperatures(tau, x, y, dt, dx, dy, cube);

  // bool intersect = true;

  // intersect= temperature_intersects_cutoff(cube);

  // return (intersect);
}
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
void SurfaceFinder::fill_cube_with_temperatures(
  Jetscape::real tau, Jetscape::real x, Jetscape::real y, 
  Jetscape::real dt, Jetscape::real dx, Jetscape::real dy,
  // std::array<std::array<std::array<double, 2>, 2>, 2>& cube
  double ***cube
  ) {
    auto tau_low = tau - dt / 2.;
    auto tau_high = tau + dt / 2.;
    auto x_left = x - dx / 2.;
    auto x_right = x + dx / 2.;
    auto y_left = y - dy / 2.;
    auto y_right = y + dy / 2.;
    
    cube[0][0][0] = bulk_info.get(tau_low, x_left, y_left, 0.0).temperature;;
    cube[0][0][1] = bulk_info.get(tau_low, x_left, y_right, 0.0).temperature;
    cube[0][1][0] = bulk_info.get(tau_low, x_right, y_left, 0.0).temperature;
    cube[0][1][1] = bulk_info.get(tau_low, x_right, y_right, 0.0).temperature;
    cube[1][0][0] = bulk_info.get(tau_high, x_left, y_left, 0.0).temperature;
    cube[1][0][1] = bulk_info.get(tau_high, x_left, y_right, 0.0).temperature;
    cube[1][1][0] = bulk_info.get(tau_high, x_right, y_left, 0.0).temperature;;
    cube[1][1][1]= bulk_info.get(tau_high, x_right, y_right, 0.0).temperature;
}
/**
 * @brief Checks if the temperature values in the cube intersect the cutoff temperature.
 * 
 * @param cube Temperature values at the corners of the cube.
 * @return True if the temperature values intersect the cutoff temperature, false otherwise.
 */
bool SurfaceFinder::temperature_intersects_cutoff(
  // const 
double ***cube
// std::array<std::array<std::array<double, 2>, 2>, 2>& cube
) {
    bool intersect = true;
    if (
        ((T_cut - cube[0][0][0]) * (cube[1][1][1] - T_cut) < 0.0) && 
        ((T_cut - cube[0][1][0]) * (cube[1][0][1] - T_cut) < 0.0) && 
        ((T_cut - cube[0][1][1]) * (cube[1][0][0] - T_cut) < 0.0) &&
        ((T_cut - cube[0][0][1]) * (cube[1][1][0] - T_cut) < 0.0)
      )
      intersect = false;

    return intersect;
}
#pragma endregion  check_intersect_3D
#pragma region Find_fill_hypersurface_3D
// a method that converting cube from std::array<std::array<std::array<double, 2>, 2>, 2> to double***
void convert_cube(std::array<std::array<std::array<double, 2>, 2>, 2>& cube, double ***cube_temp) {
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
          cube_temp[i][j][k] = cube[i][j][k];
      }
    }
  }
}

/*
* @brief Delete the local 3D array.
* a method that takes cube_temp from previous method as input and deletes it
*/
void delete_cube(double ***cube_temp) {
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++)
      delete[] cube_temp[i][j];
    delete[] cube_temp[i];
  }
  delete[] cube_temp;
}
void create_cube(double ***&cube, double init_value) {
      cube = new double **[2];
      for (int i = 0; i < 2; i++) {
        cube[i] = new double *[2];
        for (int j = 0; j < 2; j++) {
          cube[i][j] = new double[2];
          for (int k = 0; k < 2; k++)
            cube[i][j][k] = init_value;
        }
      }
    }
void SurfaceFinder::Find_full_hypersurface_3D() {
  //log that the hypersurface 3D is being runnin
  JSINFO << "Finding a 2+1D hyper-surface at T = " << T_cut << " GeV ...";
  auto grid_tau0 = bulk_info.Tau0();
  auto grid_tauf = bulk_info.TauMax();
  auto grid_x0 = bulk_info.XMin();
  auto grid_y0 = bulk_info.YMin();

  Jetscape::real grid_dt = 0.1;
  Jetscape::real grid_dx = 0.2;
  Jetscape::real grid_dy = 0.2;

  const int dim = 3;
  std::array<double, dim> lattice_spacing = {grid_dt, grid_dx, grid_dy};

  const int ntime = static_cast<int>((grid_tauf - grid_tau0) / grid_dt);
  const int nx = static_cast<int>(std::abs(2. * grid_x0) / grid_dx);
  const int ny = static_cast<int>(std::abs(2. * grid_y0) / grid_dy);

  int surface_cell_list_sz = ntime * nx * ny;
  std::vector<std::vector<SurfaceCellInfo>> surface_cell_list_local;
  surface_cell_list_local.resize(surface_cell_list_sz);

//log that we are entering the parallel region
JSINFO << "Getting into parrallel region\n";
#pragma omp parallel
{
    // std::array<std::array<std::array<double, 2>, 2>, 2> cube = {{{0.0}}};
    double ***cube;
    create_cube(cube,0.0);
    std::unique_ptr<Cornelius> cornelius_ptr(new Cornelius());
    cornelius_ptr->init(dim, T_cut, lattice_spacing.data());

#pragma omp for collapse(3)
    for (int itime = 0; itime < ntime; itime++) {
      for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
          //log that we are in the loop over time evolution
          JSINFO << "Loop over time evolution\n";

          // loop over time evolution
          auto tau_local = grid_tau0 + (itime + 0.5) * grid_dt;
          // loops over the transverse plane
          auto x_local = grid_x0 + (i + 0.5) * grid_dx;
          auto y_local = grid_y0 + (j + 0.5) * grid_dy;
          //log that we are goig to check intersect
          JSINFO << "Checking intersect\n";
          bool intersect = check_intersect_3D(tau_local, x_local, y_local,
                                              grid_dt, grid_dx, grid_dy, cube);
          //log intersect result 
          JSINFO << "Intersect result: " << intersect << "\n";
          if (intersect) {
            //log that we are in the intersect loop
            JSINFO << "Intersect loop\n";

            cornelius_ptr->find_surface_3d(cube);
            process_surface_elements(tau_local,  x_local,y_local, grid_dt, grid_dx, grid_dy, 
                                     cube, itime, nx , ny , i , j,
                                     cornelius_ptr,surface_cell_list_local);
         }
        }
      }
    }  // end of omp for loop
    delete_cube(cube);
  }  // end of parallel region
  // reduction of local 2D vector to 1D class member vector
  reduce_surface_cell_list(surface_cell_list_local, surface_cell_list_sz);
}
/*
* @brief Reduce the local 2D vector to 1D class member vector.
*
* @param surface_cell_list_local Local 2D vector.
* @param surface_cell_list_sz Size of the local 2D vector.
*/
void SurfaceFinder::reduce_surface_cell_list(std::vector<std::vector<SurfaceCellInfo>>& surface_cell_list_local, int surface_cell_list_sz){
  for (int i = 0; i < surface_cell_list_sz; i++) {
    for (int j = 0; j < surface_cell_list_local[i].size(); j++) {
      surface_cell_list.push_back(surface_cell_list_local[i][j]);
    }
  }
}
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
void SurfaceFinder::process_surface_elements(Jetscape::real tau_local, Jetscape::real x_local, Jetscape::real y_local, 
                                             Jetscape::real grid_dt, Jetscape::real grid_dx, Jetscape::real grid_dy, 
                                            //  std::array<std::array<std::array<double, 2>, 2>, 2>& cube, 
                                              double ***cube,
                                             const int itime, const int nx ,const int ny ,int i ,int j,
                                             const std::unique_ptr<Cornelius>& cornelius_ptr, 
                                             std::vector<std::vector<SurfaceCellInfo>>& surface_cell_list_local){
  // converting cube from std::array<std::array<std::array<double, 2>, 2>, 2> to double***
  //log that we are converting cube from std::array<std::array<std::array<double, 2>, 2>, 2> to double***
  // JSINFO << "Converting cube from std::array<std::array<std::array<double, 2>, 2>, 2> to double***\n";

  // double ***cube_temp = new double **[2];
  // convert_cube(cube, cube_temp) ;
  //log that we are finding the surface
  // JSINFO << "Finding the surface\n";
  // cornelius_ptr->find_surface_3d(
  //   // cube_temp
  //   cube
  //   );
  //log that we are deleting the cube
  // JSINFO << "Deleting the cube\n";
  // delete_cube(
  //   // cube_temp
  //   cube
  //   );

  for (int isurf = 0; isurf < cornelius_ptr->get_Nelements(); isurf++) {
    auto tau_center = (cornelius_ptr->get_centroid_elem(isurf, 0) +
                        tau_local - grid_dt / 2.);
    auto x_center = (cornelius_ptr->get_centroid_elem(isurf, 1) +
                      x_local - grid_dx / 2.);
    auto y_center = (cornelius_ptr->get_centroid_elem(isurf, 2) +
                      y_local - grid_dy / 2.);

    auto da_tau = cornelius_ptr->get_normal_elem(isurf, 0);
    auto da_x = cornelius_ptr->get_normal_elem(isurf, 1);
    auto da_y = cornelius_ptr->get_normal_elem(isurf, 2);

    auto fluid_cell =
        bulk_info.get(tau_center, x_center, y_center, 0.0);
    auto surface_cell =
        PrepareASurfaceCell(tau_center, x_center, y_center, 0.0,
                            da_tau, da_x, da_y, 0.0, fluid_cell);
    surface_cell_list_local[itime * nx * ny + i * ny + j].push_back(
        surface_cell);
  }                                        
}
#pragma endregion Find_fill_hypersurface_3D
bool SurfaceFinder::check_intersect_4D(Jetscape::real tau, Jetscape::real x,
                                       Jetscape::real y, Jetscape::real eta,
                                       Jetscape::real dt, Jetscape::real dx,
                                       Jetscape::real dy, Jetscape::real deta,
                                       double ****cube) {
  bool intersect = true;

  auto tau_low = tau - dt / 2.;
  auto tau_high = tau + dt / 2.;
  auto x_left = x - dx / 2.;
  auto x_right = x + dx / 2.;
  auto y_left = y - dy / 2.;
  auto y_right = y + dy / 2.;
  auto eta_left = eta - deta / 2.;
  auto eta_right = eta + deta / 2.;

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

  if ((T_cut - cube[0][0][0][0]) * (cube[1][1][1][1] - T_cut) < 0.0)
    if ((T_cut - cube[0][0][1][1]) * (cube[1][1][0][0] - T_cut) < 0.0)
      if ((T_cut - cube[0][1][0][1]) * (cube[1][0][1][0] - T_cut) < 0.0)
        if ((T_cut - cube[0][1][1][0]) * (cube[1][0][0][1] - T_cut) < 0.0)
          if ((T_cut - cube[0][0][0][1]) * (cube[1][1][1][0] - T_cut) < 0.0)
            if ((T_cut - cube[0][0][1][0]) * (cube[1][1][0][1] - T_cut) < 0.0)
              if ((T_cut - cube[0][1][0][0]) * (cube[1][0][1][1] - T_cut) < 0.0)
                if ((T_cut - cube[0][1][1][1]) * (cube[1][0][0][0] - T_cut) <
                    0.0)
                  intersect = false;

  return (intersect);
}

void SurfaceFinder::Find_full_hypersurface_4D() {
  auto grid_tau0 = bulk_info.Tau0();
  auto grid_tauf = bulk_info.TauMax();
  auto grid_x0 = bulk_info.XMin();
  auto grid_y0 = bulk_info.YMin();
  ;
  auto grid_eta0 = bulk_info.EtaMin();
  ;

  Jetscape::real grid_dt = 0.1;
  Jetscape::real grid_dx = 0.2;
  Jetscape::real grid_dy = 0.2;
  Jetscape::real grid_deta = 0.2;

  const int dim = 4;
  double lattice_spacing[dim];
  lattice_spacing[0] = grid_dt;
  lattice_spacing[1] = grid_dx;
  lattice_spacing[2] = grid_dy;
  lattice_spacing[3] = grid_deta;

  const int ntime = static_cast<int>((grid_tauf - grid_tau0) / grid_dt);
  const int nx = static_cast<int>(std::abs(2. * grid_x0) / grid_dx);
  const int ny = static_cast<int>(std::abs(2. * grid_y0) / grid_dy);
  const int neta = static_cast<int>(std::abs(2. * grid_eta0) / grid_deta);

  int surface_cell_list_sz = ntime * nx * ny * neta;
  std::vector<std::vector<SurfaceCellInfo>> surface_cell_list_local;
  surface_cell_list_local.resize(surface_cell_list_sz);

#pragma omp parallel
  {
    double ****cube = new double ***[2];
    for (int i = 0; i < 2; i++) {
      cube[i] = new double **[2];
      for (int j = 0; j < 2; j++) {
        cube[i][j] = new double *[2];
        for (int k = 0; k < 2; k++) {
          cube[i][j][k] = new double[2];
          for (int l = 0; l < 2; l++) {
            cube[i][j][k][l] = 0.0;
          }
        }
      }
    }

    std::unique_ptr<Cornelius> cornelius_ptr(new Cornelius());
    cornelius_ptr->init(dim, T_cut, lattice_spacing);

#pragma omp for collapse(4)
    for (int itime = 0; itime < ntime; itime++) {
      for (int l = 0; l < neta; l++) {
        for (int i = 0; i < nx; i++) {
          for (int j = 0; j < ny; j++) {
            // loop over time evolution
            auto tau_local = grid_tau0 + (itime + 0.5) * grid_dt;
            auto eta_local = grid_eta0 + (l + 0.5) * grid_deta;
            // loops over the transverse plane
            auto x_local = grid_x0 + (i + 0.5) * grid_dx;
            auto y_local = grid_y0 + (j + 0.5) * grid_dy;
            bool intersect =
                check_intersect_4D(tau_local, x_local, y_local, eta_local,
                                   grid_dt, grid_dx, grid_dy, grid_deta, cube);
            if (intersect) {
              cornelius_ptr->find_surface_4d(cube);
              for (int isurf = 0; isurf < cornelius_ptr->get_Nelements();
                   isurf++) {
                auto tau_center = (cornelius_ptr->get_centroid_elem(isurf, 0) +
                                   tau_local - grid_dt / 2.);
                auto x_center = (cornelius_ptr->get_centroid_elem(isurf, 1) +
                                 x_local - grid_dx / 2.);
                auto y_center = (cornelius_ptr->get_centroid_elem(isurf, 2) +
                                 y_local - grid_dy / 2.);
                auto eta_center = (cornelius_ptr->get_centroid_elem(isurf, 3) +
                                   eta_local - grid_deta / 2.);

                auto da_tau = (cornelius_ptr->get_normal_elem(isurf, 0));
                auto da_x = (cornelius_ptr->get_normal_elem(isurf, 1));
                auto da_y = (cornelius_ptr->get_normal_elem(isurf, 2));
                auto da_eta = (cornelius_ptr->get_normal_elem(isurf, 3));

                auto fluid_cell =
                    bulk_info.get(tau_center, x_center, y_center, eta_center);
                auto surface_cell = PrepareASurfaceCell(
                    tau_center, x_center, y_center, eta_center, da_tau, da_x,
                    da_y, da_eta, fluid_cell);
                surface_cell_list_local[itime * nx * ny * neta + l * nx * ny +
                                        i * ny + j]
                    .push_back(surface_cell);
              }
            }
          }
        }
      }
    }  // end of omp for loop

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        for (int k = 0; k < 2; k++) {
          delete[] cube[i][j][k];
        }
        delete[] cube[i][j];
      }
      delete[] cube[i];
    }
    delete[] cube;
  }  // end of parallel region

  // reduction of local 2D vector to 1D class member vector
  for (int i = 0; i < surface_cell_list_sz; i++) {
    for (int j = 0; j < surface_cell_list_local[i].size(); j++) {
      surface_cell_list.push_back(surface_cell_list_local[i][j]);
    }
  }
}

SurfaceCellInfo SurfaceFinder::PrepareASurfaceCell(
    Jetscape::real tau, Jetscape::real x, Jetscape::real y, Jetscape::real eta,
    Jetscape::real da0, Jetscape::real da1, Jetscape::real da2,
    Jetscape::real da3, const FluidCellInfo fluid_cell) {
  SurfaceCellInfo temp_cell;
  temp_cell.tau = tau;
  temp_cell.x = x;
  temp_cell.y = y;
  temp_cell.eta = eta;
  temp_cell.d3sigma_mu[0] = da0;
  temp_cell.d3sigma_mu[1] = da1;
  temp_cell.d3sigma_mu[2] = da2;
  temp_cell.d3sigma_mu[3] = da3;

  temp_cell.energy_density = fluid_cell.energy_density;
  temp_cell.entropy_density = fluid_cell.entropy_density;
  temp_cell.temperature = fluid_cell.temperature;
  temp_cell.pressure = fluid_cell.pressure;
  temp_cell.qgp_fraction = fluid_cell.qgp_fraction;
  temp_cell.mu_B = fluid_cell.mu_B;
  temp_cell.mu_Q = fluid_cell.mu_C;
  temp_cell.mu_S = fluid_cell.mu_S;

  double u0 =
      sqrt(1. + fluid_cell.vx * fluid_cell.vx + fluid_cell.vy * fluid_cell.vy +
           fluid_cell.vz * fluid_cell.vz);
  double uz = u0 * fluid_cell.vz;
  temp_cell.umu[0] = u0 * cosh(eta) - uz * sinh(eta);
  temp_cell.umu[1] = u0 * fluid_cell.vx;
  temp_cell.umu[2] = u0 * fluid_cell.vy;
  temp_cell.umu[3] = -u0 * sinh(eta) + uz * cosh(eta);

  temp_cell.pi[0] = fluid_cell.pi[0][0];
  temp_cell.pi[1] = fluid_cell.pi[0][1];
  temp_cell.pi[2] = fluid_cell.pi[0][2];
  temp_cell.pi[3] = fluid_cell.pi[0][3];
  temp_cell.pi[4] = fluid_cell.pi[1][1];
  temp_cell.pi[5] = fluid_cell.pi[1][2];
  temp_cell.pi[6] = fluid_cell.pi[1][3];
  temp_cell.pi[7] = fluid_cell.pi[2][2];
  temp_cell.pi[8] = fluid_cell.pi[2][3];
  temp_cell.pi[9] = fluid_cell.pi[3][3];

  temp_cell.bulk_Pi = fluid_cell.bulk_Pi;

  return (temp_cell);
}

// function that takes a vector of SurfaceCellInfo and writes it to a file
void SurfaceFinder::WriteSurfaceToFile(
    const std::vector<SurfaceCellInfo> &surface_cell_list,
    std::string filename) {
  std::ofstream file;
  file.open(filename, std::ios::app);

  for (int i = 0; i < surface_cell_list.size(); i++) {
    file << surface_cell_list[i].sfi_to_string();
  }

  file.close();
}

}  // namespace Jetscape
