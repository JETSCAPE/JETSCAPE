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

#include "JetClass.h"

#include <assert.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>

#include "JetScapeConstants.h"
#include "JetScapeLogger.h"

namespace Jetscape {

Vertex::~Vertex() { VERBOSESHOWER(9); }

}  // namespace Jetscape
