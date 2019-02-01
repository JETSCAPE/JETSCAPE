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
// -----------------------------------------
// This is a causal liquefier with the JETSCAPE framework
// -----------------------------------------

#ifndef CAUSAL_LIQUEFIER_H
#define CAUSAL_LIQUEFIER_H

#include "Liquefier.h"

class Causal_Liquefier: public Jetscape::Liquefier {
 private:

 public:
     Causal_Liquefier();
     ~Causal_Liquefier() {};

     void Init();
     void Exec();
     void Clear();

};

#endif  // CAUSAL_LIQUEFIER_H
