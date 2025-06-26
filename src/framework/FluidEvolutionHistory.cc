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
// This is a general basic class for hydrodynamics

#include "FluidEvolutionHistory.h"

#include <MakeUniqueHelper.h>
#include <string>

#include "FluidCellInfo.h"
#include "JetScapeLogger.h"
#include "LinearInterpolation.h"

namespace Jetscape {

// convert the string type entry name to enum type EntryNames
EntryName ResolveEntryName(std::string input) {
  static const std::map<std::string, EntryName> optionStrings = {
      {"energy_density", ENTRY_ENERGY_DENSITY},
      {"entropy_density", ENTRY_ENTROPY_DENSITY},
      {"temperature", ENTRY_TEMPERATURE},
      {"pressure", ENTRY_PRESSURE},
      {"qgp_fraction", ENTRY_QGP_FRACTION},
      {"mu_b", ENTRY_MU_B},
      {"mu_c", ENTRY_MU_C},
      {"mu_s", ENTRY_MU_S},
      {"vx", ENTRY_VX},
      {"vy", ENTRY_VY},
      {"vz", ENTRY_VZ},
      {"pi00", ENTRY_PI00},
      {"pi01", ENTRY_PI01},
      {"pi02", ENTRY_PI02},
      {"pi03", ENTRY_PI03},
      {"pi11", ENTRY_PI11},
      {"pi12", ENTRY_PI12},
      {"pi13", ENTRY_PI13},
      {"pi22", ENTRY_PI22},
      {"pi23", ENTRY_PI23},
      {"pi33", ENTRY_PI33},
      {"bulk_pi", ENTRY_BULK_PI},
  };
  auto itr = optionStrings.find(input);
  if (itr != optionStrings.end()) {
    return itr->second;
  } else {
    return ENTRY_INVALID;
  }
}

// It checks whether a space-time point (tau, x, y, eta) is inside evolution
// history or outside.
int EvolutionHistory::CheckInRange(Jetscape::real tau, Jetscape::real x,
                                   Jetscape::real y, Jetscape::real eta) const {
  int status = 1;
  if (tau < tau_min || tau > TauMax()) {
    std::string warn_message =
        ("tau=" + std::to_string(tau) + " is not in range [" +
         std::to_string(tau_min) + "," + std::to_string(TauMax()) + "]");
    // throw InvalidSpaceTimeRange(warn_message);
    // JSWARN << warn_message;
    status = 0;
  }
  if (x < x_min || x > XMax()) {
    std::string warn_message =
        ("x=" + std::to_string(x) + " is not in range [" +
         std::to_string(x_min) + "," + std::to_string(XMax()) + "]");
    // throw InvalidSpaceTimeRange(warn_message);
    // JSWARN << warn_message;
    status = 0;
  }
  if (y < y_min || y > YMax()) {
    std::string warn_message =
        ("y=" + std::to_string(y) + " is not in range [" +
         std::to_string(y_min) + "," + std::to_string(YMax()) + "]");
    // throw InvalidSpaceTimeRange(warn_message);
    // JSWARN << warn_message;
    status = 0;
  }
  if (!boost_invariant) {
    if (eta < eta_min || eta > EtaMax()) {
      std::string warn_message =
          ("eta=" + std::to_string(eta) + " is not in range [" +
           std::to_string(eta_min) + "," + std::to_string(EtaMax()) + "]");
      // throw InvalidSpaceTimeRange(warn_message);
      // JSWARN << warn_message;
      status = 0;
    }
  }
  return (status);
}

/** Construct evolution history given the bulk_data and the data_info */
void EvolutionHistory::FromVector(const std::vector<float> &data_,
                                  const std::vector<std::string> &data_info_,
                                  float tau_min_, float dtau_, float x_min_,
                                  float dx_, int nx_, float y_min_, float dy_,
                                  int ny_, float eta_min_, float deta_,
                                  int neta_, bool tau_eta_is_tz_) {
  data_vector = data_;
  data_info = data_info_;
  tau_min = tau_min_;
  x_min = x_min_;
  y_min = y_min_;
  eta_min = eta_min_;
  dtau = dtau_;
  dx = dx_;
  dy = dy_;
  deta = deta_;
  nx = nx_;
  ny = ny_;
  neta = neta_;
  tau_eta_is_tz = tau_eta_is_tz_;
  ntau = data_.size() / (data_info_.size() * nx * ny * neta);
}

/* This function will read the sparse data stored in data_ with associated
 * information data_info_ into to FluidCellInfo object */
FluidCellInfo EvolutionHistory::GetFluidCell(int id_tau, int id_x, int id_y,
                                             int id_eta) const {
  int entries_per_record = data_info.size();
  int id_eta_corrected = id_eta;
  // set id_eta=0 if hydro is in 2+1D mode
  if (neta == 0 || neta == 1) {
    id_eta_corrected = 0;
  }

  int record_starting_id = CellIndex(id_tau, id_x, id_y, id_eta_corrected);

  // if data_vector and data_info are not used to construct evolution history
  // then the data should have the format of vector<FluidCellInfo>.
  if (entries_per_record == 0) {
    return data.at(record_starting_id);
  }
  // otherwise construct the fluid cell info from data_vector and data_info
  auto fluid_cell_ptr = make_unique<FluidCellInfo>();

  record_starting_id *= entries_per_record;
  for (int i = 0; i < entries_per_record; i++) {
    auto entry_name = ResolveEntryName(data_info.at(i));
    auto entry_data = data_vector.at(record_starting_id + i);
    switch (entry_name) {
      case ENTRY_ENERGY_DENSITY:
        fluid_cell_ptr->energy_density = entry_data;
        break;
      case ENTRY_ENTROPY_DENSITY:
        fluid_cell_ptr->entropy_density = entry_data;
        break;
      case ENTRY_TEMPERATURE:
        fluid_cell_ptr->temperature = entry_data;
        break;
      case ENTRY_PRESSURE:
        fluid_cell_ptr->pressure = entry_data;
        break;
      case ENTRY_QGP_FRACTION:
        fluid_cell_ptr->qgp_fraction = entry_data;
        break;
      case ENTRY_MU_B:
        fluid_cell_ptr->mu_B = entry_data;
        break;
      case ENTRY_MU_C:
        fluid_cell_ptr->mu_C = entry_data;
        break;
      case ENTRY_MU_S:
        fluid_cell_ptr->mu_S = entry_data;
        break;
      case ENTRY_VX:
        fluid_cell_ptr->vx = entry_data;
        break;
      case ENTRY_VY:
        fluid_cell_ptr->vy = entry_data;
        break;
      case ENTRY_VZ:
        fluid_cell_ptr->vz = entry_data;
        break;
      case ENTRY_PI00:
        fluid_cell_ptr->pi[0][0] = entry_data;
        break;
      case ENTRY_PI01:
        fluid_cell_ptr->pi[0][1] = entry_data;
        fluid_cell_ptr->pi[1][0] = entry_data;
        break;
      case ENTRY_PI02:
        fluid_cell_ptr->pi[0][2] = entry_data;
        fluid_cell_ptr->pi[2][0] = entry_data;
        break;
      case ENTRY_PI03:
        fluid_cell_ptr->pi[0][3] = entry_data;
        fluid_cell_ptr->pi[3][0] = entry_data;
        break;
      case ENTRY_PI11:
        fluid_cell_ptr->pi[1][1] = entry_data;
        break;
      case ENTRY_PI12:
        fluid_cell_ptr->pi[1][2] = entry_data;
        fluid_cell_ptr->pi[2][1] = entry_data;
        break;
      case ENTRY_PI13:
        fluid_cell_ptr->pi[1][3] = entry_data;
        fluid_cell_ptr->pi[3][1] = entry_data;
        break;
      case ENTRY_PI22:
        fluid_cell_ptr->pi[2][2] = entry_data;
        break;
      case ENTRY_PI23:
        fluid_cell_ptr->pi[2][3] = entry_data;
        fluid_cell_ptr->pi[3][2] = entry_data;
        break;
      case ENTRY_PI33:
        fluid_cell_ptr->pi[3][3] = entry_data;
        break;
      case ENTRY_BULK_PI:
        fluid_cell_ptr->bulk_Pi = entry_data;
        break;
      default:
        JSWARN << "The entry name in data_info_ must be one of the \
                        energy_density, entropy_density, temperature, pressure, qgp_fraction, \
                        mu_b, mu_c, mu_s, vx, vy, vz, pi00, pi01, pi02, pi03, pi11, pi12, \
                        pi13, pi22, pi23, pi33, bulk_pi";
        break;
    }
  }

  return *fluid_cell_ptr;
}

/** For one given time step id_tau,
 * get FluidCellInfo at spatial point (x, y, eta)*/
FluidCellInfo EvolutionHistory::GetAtTimeStep(int id_tau, Jetscape::real x,
                                              Jetscape::real y,
                                              Jetscape::real eta) const {
  int id_x = GetIdX(x);
  int id_y = GetIdY(y);
  int id_eta = 0;
  if (!boost_invariant)
    id_eta = GetIdEta(eta);

  auto c000 = GetFluidCell(id_tau, id_x, id_y, id_eta);
  auto c001 = GetFluidCell(id_tau, id_x, id_y, id_eta + 1);
  auto c010 = GetFluidCell(id_tau, id_x, id_y + 1, id_eta);
  auto c011 = GetFluidCell(id_tau, id_x, id_y + 1, id_eta + 1);
  auto c100 = GetFluidCell(id_tau, id_x + 1, id_y, id_eta);
  auto c101 = GetFluidCell(id_tau, id_x + 1, id_y, id_eta + 1);
  auto c110 = GetFluidCell(id_tau, id_x + 1, id_y + 1, id_eta);
  auto c111 = GetFluidCell(id_tau, id_x + 1, id_y + 1, id_eta + 1);
  real x0 = XCoord(id_x);
  real x1 = XCoord(id_x + 1);
  real y0 = YCoord(id_y);
  real y1 = YCoord(id_y + 1);
  real eta0 = EtaCoord(id_eta);
  auto eta1 = 0.0;
  if (!boost_invariant)
    eta1 = EtaCoord(id_eta + 1);

  return TrilinearInt(x0, x1, y0, y1, eta0, eta1, c000, c001, c010, c011, c100,
                      c101, c110, c111, x, y, eta);
}

// do interpolation along time direction; we may also need high order
// interpolation functions
FluidCellInfo EvolutionHistory::get(Jetscape::real tau, Jetscape::real x,
                                    Jetscape::real y,
                                    Jetscape::real eta) const {
  int status = CheckInRange(tau, x, y, eta);
  if (status == 0) {
    FluidCellInfo zero_cell;
    return (zero_cell);
  }
  int id_tau = GetIdTau(tau);
  auto tau0 = TauCoord(id_tau);
  auto tau1 = TauCoord(id_tau + 1);
  auto bulk0 = GetAtTimeStep(id_tau, x, y, eta);
  auto bulk1 = GetAtTimeStep(id_tau + 1, x, y, eta);
  return (LinearInt(tau0, tau1, bulk0, bulk1, tau));
}

FluidCellInfo EvolutionHistory::get_tz(Jetscape::real t, Jetscape::real x,
                                       Jetscape::real y,
                                       Jetscape::real z) const {
  Jetscape::real tau = 0.0;
  Jetscape::real eta = 0.0;
  if (t * t > z * z) {
    tau = sqrt(t * t - z * z);
    eta = 0.5 * log((t + z) / (t - z));
  } else {
    VERBOSE(4) << "the quest point is outside the light cone! "
               << "t = " << t << ", z = " << z;
  }
  auto cell = get(tau, x, y, eta);
  if (boost_invariant) {
    cell.vz = z / t;
    Jetscape::real gammaL = 1.0 / sqrt(1.0 - cell.vz * cell.vz);
    cell.vx /= gammaL;
    cell.vy /= gammaL;
  }
  return (cell);
}

}  // end namespace Jetscape
