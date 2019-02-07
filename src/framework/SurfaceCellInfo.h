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

#ifndef SURFACECELLINFO_H
#define SURFACECELLINFO_H

#include "RealType.h"

namespace Jetscape {

class SurfaceCellInfo {
 public:
    // data structure for outputing hyper-surface information
    Jetscape::real tau;
    Jetscape::real x;
    Jetscape::real y;
    Jetscape::real eta;
    Jetscape::real d3sigma_mu[4];     //!< Surface vector.
    Jetscape::real energy_density;    //!< Local energy density [GeV/fm^3].
    Jetscape::real entropy_density;   //!< Local entropy density [1/fm^3].
    Jetscape::real temperature;       //!< Local temperature [GeV].
    Jetscape::real pressure;          //!< Thermal pressure [GeV/fm^3].
    Jetscape::real qgp_fraction;      //!< Fraction of quark gluon plasma assuming medium is in QGP+HRG phase.
    Jetscape::real mu_B;              //!< Net baryon chemical potential [GeV].
    Jetscape::real mu_C;              //!< Net charge chemical potential [GeV].
    Jetscape::real mu_S;              //!< Net strangeness chemical potential [GeV].
    Jetscape::real vx, vy, vz;        //!< Flow velocity.
    Jetscape::real pi[4][4];          //!< Shear stress tensor [GeV/fm^3].
    Jetscape::real bulk_Pi;           //!< Bulk viscous pressure [GeV/fm^3].
    
    /** Default constructor. */
    SurfaceCellInfo() = default;
    
    /** Destructor. */
    ~SurfaceCellInfo() {};
};

}

#endif  // SURFACECELLINFO_H
