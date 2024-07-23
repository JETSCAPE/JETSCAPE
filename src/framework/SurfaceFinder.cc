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
#include <array>
#include "RealType.h"
#include "SurfaceFinder.h"
#include "cornelius.h"
#include "FluidEvolutionHistory.h"
#include "JetScapeLogger.h"

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
  if (boost_invariant) {
    JSINFO << "Finding a 2+1D hyper-surface at T = " << T_cut << " GeV ...";
    Find_full_hypersurface_3D();
  } else {
    JSINFO << "Finding a 3+1D hyper-surface at T = " << T_cut << " GeV ...";
    Find_full_hypersurface_4D();
  }
}

bool SurfaceFinder::check_intersect_3D(Jetscape::real tau, Jetscape::real x,
                                       Jetscape::real y, Jetscape::real dt,
                                       Jetscape::real dx, Jetscape::real dy,
                                       double ***cube) {
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
}

void SurfaceFinder::Find_full_hypersurface_3D() {
  auto grid_tau0 = bulk_info.Tau0();
  auto grid_tauf = bulk_info.TauMax();
  auto grid_x0 = bulk_info.XMin();
  auto grid_y0 = bulk_info.YMin();
  ;

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

  const int ntime = static_cast<int>((grid_tauf - grid_tau0) / grid_dt);
  const int nx = static_cast<int>(std::abs(2. * grid_x0) / grid_dx);
  const int ny = static_cast<int>(std::abs(2. * grid_y0) / grid_dy);

  double ***cube = new double **[2];
  for (int i = 0; i < 2; i++) {
    cube[i] = new double *[2];
    for (int j = 0; j < 2; j++) {
      cube[i][j] = new double[2];
      for (int k = 0; k < 2; k++)
        cube[i][j][k] = 0.0;
    }
  }

  for (int itime = 0; itime < ntime; itime++) {
    // loop over time evolution
    auto tau_local = grid_tau0 + (itime + 0.5) * grid_dt;
    for (int i = 0; i < nx; i++) {
      // loops over the transverse plane
      auto x_local = grid_x0 + (i + 0.5) * grid_dx;
      for (int j = 0; j < ny; j++) {
        auto y_local = grid_y0 + (j + 0.5) * grid_dy;
        bool intersect = check_intersect_3D(tau_local, x_local, y_local,
                                            grid_dt, grid_dx, grid_dy, cube);
        if (intersect) {
          cornelius_ptr->find_surface_3d(cube);
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
                PrepareASurfaceCell(tau_center, x_center, y_center, 0.0, da_tau,
                                    da_x, da_y, 0.0, fluid_cell);
            surface_cell_list.push_back(surface_cell);
          }
        }
      }
    }
  }

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++)
      delete[] cube[i][j];
    delete[] cube[i];
  }
  delete[] cube;
}


#pragma region check intersect 4D
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
bool SurfaceFinder::check_intersect_4D(
    Jetscape::real tau, Jetscape::real x, Jetscape::real y, Jetscape::real eta,
    Jetscape::real dt, Jetscape::real dx, Jetscape::real dy, Jetscape::real deta,
    std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> &cube) {
    
    
    fill_cube_with_temperatures(tau, x, y, eta, dt, dx, dy, deta, cube);
    
    bool intersects=true;
    intersects=!temperature_intersects_cutoff(cube);

    return intersects;
}

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
void SurfaceFinder::fill_cube_with_temperatures(
    Jetscape::real tau, Jetscape::real x, Jetscape::real y, Jetscape::real eta,
    Jetscape::real dt, Jetscape::real dx, Jetscape::real dy, Jetscape::real deta,
    std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> &cube) {

    auto tau_low = tau - dt / 2.;
    auto tau_high = tau + dt / 2.;
    auto x_left = x - dx / 2.;
    auto x_right = x + dx / 2.;
    auto y_left = y - dy / 2.;
    auto y_right = y + dy / 2.;
    auto eta_left = eta - deta / 2.;
    auto eta_right = eta + deta / 2.;

    cube[0][0][0][0] = bulk_info.get(tau_low, x_left, y_left, eta_left).temperature;
    cube[0][0][0][1] = bulk_info.get(tau_low, x_left, y_left, eta_right).temperature;
    cube[0][0][1][0] = bulk_info.get(tau_low, x_left, y_right, eta_left).temperature;
    cube[0][0][1][1] = bulk_info.get(tau_low, x_left, y_right, eta_right).temperature;
    cube[0][1][0][0] = bulk_info.get(tau_low, x_right, y_left, eta_left).temperature;
    cube[0][1][0][1] = bulk_info.get(tau_low, x_right, y_left, eta_right).temperature;
    cube[0][1][1][0] = bulk_info.get(tau_low, x_right, y_right, eta_left).temperature;
    cube[0][1][1][1] = bulk_info.get(tau_low, x_right, y_right, eta_right).temperature;
    cube[1][0][0][0] = bulk_info.get(tau_high, x_left, y_left, eta_left).temperature;
    cube[1][0][0][1] = bulk_info.get(tau_high, x_left, y_left, eta_right).temperature;
    cube[1][0][1][0] = bulk_info.get(tau_high, x_left, y_right, eta_left).temperature;
    cube[1][0][1][1] = bulk_info.get(tau_high, x_left, y_right, eta_right).temperature;
    cube[1][1][0][0] = bulk_info.get(tau_high, x_right, y_left, eta_left).temperature;
    cube[1][1][0][1] = bulk_info.get(tau_high, x_right, y_left, eta_right).temperature;
    cube[1][1][1][0] = bulk_info.get(tau_high, x_right, y_right, eta_left).temperature;
    cube[1][1][1][1] = bulk_info.get(tau_high, x_right, y_right, eta_right).temperature;
}

/**
 * @brief Checks if the temperature in the 4D array cube intersects the cutoff temperature.
 *
 * @param cube 4D array containing temperature values of the grid cell.
 * @return True if the temperature intersects the cutoff, false otherwise.
 */
bool SurfaceFinder::temperature_intersects_cutoff(
    const std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> &cube) {

    return (T_cut - cube[0][0][0][0]) * (cube[1][1][1][1] - T_cut) >= 0.0 &&
           (T_cut - cube[0][0][1][1]) * (cube[1][1][0][0] - T_cut) >= 0.0 &&
           (T_cut - cube[0][1][0][1]) * (cube[1][0][1][0] - T_cut) >= 0.0 &&
           (T_cut - cube[0][1][1][0]) * (cube[1][0][0][1] - T_cut) >= 0.0 &&
           (T_cut - cube[0][0][0][1]) * (cube[1][1][1][0] - T_cut) >= 0.0 &&
           (T_cut - cube[0][0][1][0]) * (cube[1][1][0][1] - T_cut) >= 0.0 &&
           (T_cut - cube[0][1][0][0]) * (cube[1][0][1][1] - T_cut) >= 0.0 &&
           (T_cut - cube[0][1][1][1]) * (cube[1][0][0][0] - T_cut) >= 0.0;
}
#pragma endregion check intersect 4D

#pragma region  finding full hypersurface 4D
void SurfaceFinder::Find_full_hypersurface_4D() {
    auto grid_tau0 = bulk_info.Tau0();
    auto grid_tauf = bulk_info.TauMax();
    auto grid_x0 = bulk_info.XMin();
    auto grid_y0 = bulk_info.YMin();
    auto grid_eta0 = bulk_info.EtaMin();

    const Jetscape::real grid_dt = 0.1;
    const Jetscape::real grid_dx = 0.2;
    const Jetscape::real grid_dy = 0.2;
    const Jetscape::real grid_deta = 0.2;

    const int dim = 4;
    const double lattice_spacing[dim] = {grid_dt, grid_dx, grid_dy, grid_deta};

    std::unique_ptr<Cornelius> cornelius_ptr(new Cornelius());
    cornelius_ptr->init(dim, T_cut, lattice_spacing);

    const int ntime = static_cast<int>((grid_tauf - grid_tau0) / grid_dt);
    const int nx = static_cast<int>(std::abs(2. * grid_x0) / grid_dx);
    const int ny = static_cast<int>(std::abs(2. * grid_y0) / grid_dy);
    const int neta = static_cast<int>(std::abs(2. * grid_eta0) / grid_deta);

    auto cube = create_cube();

    for (int itime = 0; itime < ntime; ++itime) {
        process_time_slice(itime, grid_tau0, grid_dt, grid_eta0, grid_deta, neta, grid_x0, grid_dx, nx, grid_y0, grid_dy, ny, cube, cornelius_ptr);
    }
}

std::vector<std::vector<std::vector<std::vector<double>>>> SurfaceFinder::create_cube() {
    return std::vector<std::vector<std::vector<std::vector<double>>>>(
        2, std::vector<std::vector<std::vector<double>>>(
               2, std::vector<std::vector<double>>(
                      2, std::vector<double>(2, 0.0))));
}

void SurfaceFinder::process_time_slice(
    int itime, Jetscape::real grid_tau0, Jetscape::real grid_dt,
    Jetscape::real grid_eta0, Jetscape::real grid_deta, int neta,
    Jetscape::real grid_x0, Jetscape::real grid_dx, int nx,
    Jetscape::real grid_y0, Jetscape::real grid_dy, int ny,
    std::vector<std::vector<std::vector<std::vector<double>>>>& cube,
    const std::unique_ptr<Cornelius>& cornelius_ptr) {

    const auto tau_local = grid_tau0 + (itime + 0.5) * grid_dt;
    for (int l = 0; l < neta; ++l) {
        process_eta_slice(l, tau_local, grid_eta0, grid_deta, grid_x0, grid_dx, nx, grid_y0, grid_dy, ny, cube, cornelius_ptr);
    }
}

void SurfaceFinder::process_eta_slice(
    int l, Jetscape::real tau_local, Jetscape::real grid_eta0, Jetscape::real grid_deta,
    Jetscape::real grid_x0, Jetscape::real grid_dx, int nx,
    Jetscape::real grid_y0, Jetscape::real grid_dy, int ny,
    std::vector<std::vector<std::vector<std::vector<double>>>>& cube,
    const std::unique_ptr<Cornelius>& cornelius_ptr) {

    const auto eta_local = grid_eta0 + (l + 0.5) * grid_deta;
    for (int i = 0; i < nx; ++i) {
        process_x_slice(i, tau_local, eta_local, grid_x0, grid_dx, grid_y0, grid_dy, ny, cube, cornelius_ptr);
    }
}

void SurfaceFinder::process_x_slice(
    int i, Jetscape::real tau_local, Jetscape::real eta_local, Jetscape::real grid_x0, Jetscape::real grid_dx,
    Jetscape::real grid_y0, Jetscape::real grid_dy, int ny,
    std::vector<std::vector<std::vector<std::vector<double>>>>& cube,
    const std::unique_ptr<Cornelius>& cornelius_ptr) {

    const auto x_local = grid_x0 + (i + 0.5) * grid_dx;
    for (int j = 0; j < ny; ++j) {
        const auto y_local = grid_y0 + (j + 0.5) * grid_dy;
        process_grid_point(tau_local, x_local, y_local, eta_local, grid_dt, grid_dx, grid_dy, grid_deta, cube, cornelius_ptr);
    }
}

void SurfaceFinder::process_grid_point(
    Jetscape::real tau_local, Jetscape::real x_local, Jetscape::real y_local, Jetscape::real eta_local,
    Jetscape::real grid_dt, Jetscape::real grid_dx, Jetscape::real grid_dy, Jetscape::real grid_deta,
    std::vector<std::vector<std::vector<std::vector<double>>>>& cube,
    const std::unique_ptr<Cornelius>& cornelius_ptr) {

    bool intersect = check_intersect_4D(tau_local, x_local, y_local, eta_local, grid_dt, grid_dx, grid_dy, grid_deta, cube);
    if (intersect) {
        cornelius_ptr->find_surface_4d(cube);
        for (int isurf = 0; isurf < cornelius_ptr->get_Nelements(); ++isurf) {
            auto [tau_center, x_center, y_center, eta_center] = compute_centroids(cornelius_ptr, isurf, tau_local, x_local, y_local, eta_local, grid_dt, grid_dx, grid_dy, grid_deta);
            auto [da_tau, da_x, da_y, da_eta] = compute_normals(cornelius_ptr, isurf);
            auto fluid_cell = bulk_info.get(tau_center, x_center, y_center, eta_center);
            auto surface_cell = PrepareASurfaceCell(tau_center, x_center, y_center, eta_center, da_tau, da_x, da_y, da_eta, fluid_cell);
            surface_cell_list.push_back(surface_cell);
        }
    }
}

std::tuple<Jetscape::real, Jetscape::real, Jetscape::real, Jetscape::real> SurfaceFinder::compute_centroids(
    const std::unique_ptr<Cornelius>& cornelius_ptr, int isurf, Jetscape::real tau_local, Jetscape::real x_local, Jetscape::real y_local, Jetscape::real eta_local,
    Jetscape::real grid_dt, Jetscape::real grid_dx, Jetscape::real grid_dy, Jetscape::real grid_deta) {

    auto tau_center = cornelius_ptr->get_centroid_elem(isurf, 0) + tau_local - grid_dt / 2.;
    auto x_center = cornelius_ptr->get_centroid_elem(isurf, 1) + x_local - grid_dx / 2.;
    auto y_center = cornelius_ptr->get_centroid_elem(isurf, 2) + y_local - grid_dy / 2.;
    auto eta_center = cornelius_ptr->get_centroid_elem(isurf, 3) + eta_local - grid_deta / 2.;
    return {tau_center, x_center, y_center, eta_center};
}

std::tuple<Jetscape::real, Jetscape::real, Jetscape::real, Jetscape::real> SurfaceFinder::compute_normals(
    const std::unique_ptr<Cornelius>& cornelius_ptr, int isurf) {

    auto da_tau = cornelius_ptr->get_normal_elem(isurf, 0);
    auto da_x = cornelius_ptr->get_normal_elem(isurf, 1);
    auto da_y = cornelius_ptr->get_normal_elem(isurf, 2);
    auto da_eta = cornelius_ptr->get_normal_elem(isurf, 3);
    return {da_tau, da_x, da_y, da_eta};
}

#pragma endregion




SurfaceCellInfo SurfaceFinder::PrepareASurfaceCell(
    Jetscape::real tau, Jetscape::real x, Jetscape::real y, Jetscape::real eta,
    Jetscape::real da0, Jetscape::real da1, Jetscape::real da2,
    Jetscape::real da3, const FluidCellInfo& fluid_cell) {

  SurfaceCellInfo temp_cell {
    .tau = tau,
    .x = x,
    .y = y,
    .eta = eta,
    .d3sigma_mu = {da0, da1, da2, da3},
    .energy_density = fluid_cell.energy_density,
    .entropy_density = fluid_cell.entropy_density,
    .temperature = fluid_cell.temperature,
    .pressure = fluid_cell.pressure,
    .qgp_fraction = fluid_cell.qgp_fraction,
    .mu_B = fluid_cell.mu_B,
    .mu_Q = fluid_cell.mu_Q,
    .mu_S = fluid_cell.mu_S,
  };

  const auto [vx, vy, vz] = std::tie(fluid_cell.vx, fluid_cell.vy, fluid_cell.vz);
  const auto u0 = std::sqrt(1. + vx * vx + vy * vy + vz * vz);
  const auto uz = u0 * vz;

  temp_cell.umu[0] = u0 * std::cosh(eta) - uz * std::sinh(eta);
  temp_cell.umu[1] = u0 * vx;
  temp_cell.umu[2] = u0 * vy;
  temp_cell.umu[3] = -u0 * std::sinh(eta) + uz * std::cosh(eta);

  std::copy(&fluid_cell.pi[0][0], &fluid_cell.pi[0][0] + 10, &temp_cell.pi[0]);

  temp_cell.bulk_Pi = fluid_cell.bulk_Pi;

  return temp_cell;
}


} // namespace Jetscape
