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

#ifndef JETSCAPECONSTANTS_H
#define JETSCAPECONSTANTS_H

namespace Jetscape {

/// Constant: Pi value
static double pi = 3.141592653589793;

/// Constant: Number of flavors
static double nf = 3.0;

/// Constant: Color factor associated with gluon emission from a quark
static double Cf = 4.0 / 3.0;

/// Constant: XXX
static double Tf = 0.5;

/// Constant: Color factor associated with gluon emission from a gluon
static double Ca = 3.0;

/// Constant: Number of colors
static double Nc = 3.0;

/// Constant: QCD scale in GeV, 0.4 is chosen in JETSET
static double Lambda_QCD = 0.2;

/// Constant: hbarC in GeV fm
static const double hbarC = 0.197327053;

/// Constant: Conversion factor from fm to GeV^-1
static const double fmToGeVinv = 1.0 / hbarC;

static double zeta3 = 1.20206;

static double mu = 0.722;

/// Constant:: Maximum value from the standard C++ random number generator
static double maxN = double(pow(2.0, 31.0) - 1.0);

/// Constant: Maximum double value
static double a_very_large_number = maxN;

/* the following 2 lines control the error in the analytical part of the
 * calculation. */
/* Note analytical approximation, cannot be rectified by more statistics */
/* However, more accurate analytical calculation will require less statistics to
 * obtain smooth distributions */
/* the value is something for the user to choose based on his/her computing
 * resources  */

static double error = 0.02;
static double approx = 0.02;

static double s_error = 0.01;
static double s_approx = 0.01;

static double E_minimum = 1.0;

/// Constant: Rounding error, slightly more than float precision
static double rounding_error = 1e-6;

/******************************************************************************/
/// Constant: Gluon PDG code
static int gid = 21;

/// Constant: Down quark PDG code
static int qid = 1;

/// Constant: Up quark PDG code
static int uid = 2;

/// Constant: Strange quark PDG code
static int did = 1;

/// Constant: Charm quark PDG code
static int sid = 3;

/// Constant: Bottom quark PDG code
static int cid = 4;

/// Constant: Top quark PDG code
static int bid = 5;

/// Constant: Photon PDG code
static int photonid = 22;

/// Constant: Z boson PDG code
static int Zid = 23;

/// Constant: W+ boson PDG code, use -24 for W-
static int Wid = 24;
/******************************************************************************/

};      // namespace Jetscape
#endif  // JETSCAPECONSTANTS_H
