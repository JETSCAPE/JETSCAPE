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

// This is a general configuration for all files

#ifndef REALTYPE_H
#define REALTYPE_H
#include <tuple>

namespace Jetscape {

typedef float real;
typedef std::tuple<real, real, real> real3;
typedef std::tuple<real, real, real, real> real4;

} // end namespace Jetscape

#endif  // REALTYPE_H
