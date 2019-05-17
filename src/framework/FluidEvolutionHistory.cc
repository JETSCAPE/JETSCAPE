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
// This is a general basic class for hydrodynamics

#include <string>
#include "FluidEvolutionHistory.h"
#include "FluidCellInfo.h"
#include "LinearInterpolation.h"
#include "JetScapeLogger.h"

namespace Jetscape {

// It checks whether a space-time point (tau, x, y, eta) is inside evolution
// history or outside.
int EvolutionHistory::CheckInRange(Jetscape::real tau, Jetscape::real x,
                                Jetscape::real y, Jetscape::real eta) const {
    int status = 1;
    if (tau < tau_min || tau > TauMax()) {
        std::string warn_message = ("tau=" + std::to_string(tau)
                + " is not in range [" + std::to_string(tau_min) + ","
                + std::to_string(TauMax()) + "]");
        //throw InvalidSpaceTimeRange(warn_message);
        //JSWARN << warn_message;
        status = 0;
    }
    if (x < x_min || x > XMax()) {
        std::string warn_message = ("x=" + std::to_string(x)
                + " is not in range [" + std::to_string(x_min) + ","
                + std::to_string(XMax()) + "]");
        //throw InvalidSpaceTimeRange(warn_message);
        //JSWARN << warn_message;
        status = 0;
    }
    if (y < y_min || y > YMax()) {
        std::string warn_message = ("y=" + std::to_string(y)
                + " is not in range [" + std::to_string(y_min) + "," 
                + std::to_string(YMax()) + "]");
        //throw InvalidSpaceTimeRange(warn_message);
        //JSWARN << warn_message;
        status = 0;
    }
    if (!boost_invariant) {
        if (eta < eta_min || eta > EtaMax()) {
            std::string warn_message = ("eta=" + std::to_string(eta)
                    + " is not in range [" + std::to_string(eta_min) + "," 
                    + std::to_string(EtaMax()) + "]");
            //throw InvalidSpaceTimeRange(warn_message);
            //JSWARN << warn_message;
            status = 0;
        }
    }
    return(status);
}
  
/** For one given time step id_tau,
   * get FluidCellInfo at spatial point (x, y, eta)*/
FluidCellInfo EvolutionHistory::GetAtTimeStep(
        int id_tau, Jetscape::real x, Jetscape::real y,
        Jetscape::real eta) const {
    int id_x   = GetIdX(x);
    int id_y   = GetIdY(y);
    int id_eta = 0;
    if (!boost_invariant)
        id_eta = GetIdEta(eta);

    // cijk for idx=i, idy=j and id_eta=k
    int c000 = CellIndex(id_tau, id_x, id_y, id_eta);
    int c001 = CellIndex(id_tau, id_x, id_y, id_eta+1);
    int c010 = CellIndex(id_tau, id_x, id_y+1, id_eta);
    int c011 = CellIndex(id_tau, id_x, id_y+1, id_eta+1);
    int c100 = CellIndex(id_tau, id_x+1, id_y, id_eta);
    int c101 = CellIndex(id_tau, id_x+1, id_y, id_eta+1);
    int c110 = CellIndex(id_tau, id_x+1, id_y+1, id_eta);
    int c111 = CellIndex(id_tau, id_x+1, id_y+1, id_eta+1);
    auto x0 = XCoord(id_x);
    auto x1 = XCoord(id_x + 1);
    auto y0 = YCoord(id_y);
    auto y1 = YCoord(id_y + 1);
    auto eta0 = EtaCoord(id_eta);
    auto eta1 = 0.0;
    if (!boost_invariant)
        eta1 = EtaCoord(id_eta + 1);

    return(TrilinearInt(x0, x1, y0, y1, eta0, eta1,
                data.at(c000), data.at(c001), data.at(c010), data.at(c011),
                data.at(c100), data.at(c101), data.at(c110), data.at(c111),
                x, y, eta));
}
  
// do interpolation along time direction; we may also need high order
// interpolation functions 
FluidCellInfo EvolutionHistory::get(Jetscape::real tau, Jetscape::real x,
                                    Jetscape::real y,
                                    Jetscape::real eta) const {
    int status = CheckInRange(tau, x, y, eta);
    if (status == 0) {
        FluidCellInfo zero_cell;
        return(zero_cell);
    }
    int id_tau = GetIdTau(tau);
    auto tau0  = TauCoord(id_tau);
    auto tau1  = TauCoord(id_tau + 1);
    auto bulk0 = GetAtTimeStep(id_tau, x, y, eta);
    auto bulk1 = GetAtTimeStep(id_tau + 1, x, y, eta);
    return(LinearInt(tau0, tau1, bulk0, bulk1, tau));
}
    
FluidCellInfo EvolutionHistory::get_tz(Jetscape::real t, Jetscape::real x,
                                       Jetscape::real y,
                                       Jetscape::real z) const {
    Jetscape::real tau = 0.0;
    Jetscape::real eta = 0.0;
    if (t*t > z*z) {
        tau = sqrt(t*t - z*z);
        eta = 0.5*log((t + z)/(t - z));
    } else {
        JSWARN << "the quest point is outside the light cone! "
               << "t = " << t << ", z = " << z;
    }
    return(get(tau, x, y, eta));
}

}  // end namespace Jetscape
