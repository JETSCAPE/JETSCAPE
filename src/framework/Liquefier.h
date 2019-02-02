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

#ifndef LIQUEFIER_H
#define LIQUEFIER_H

#include "JetScapeModuleBase.h"
#include "RealType.h"

namespace Jetscape {


class Liquefier : public JetScapeModuleBase {
 private:

 public:
     Liquefier();
     ~Liquefier() {};

     void Init();
     void Exec();

     void get_source(Jetscape::real t, Jetscape::real &j0);

};

};

#endif  // LIQUEFIER_H

