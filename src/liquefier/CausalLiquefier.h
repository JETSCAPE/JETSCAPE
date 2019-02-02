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

#ifndef CAUSALLIQUEFIER_H
#define CAUSALLIQUEFIER_H

#include "Liquefier.h"

class CausalLiquefier: public Jetscape::Liquefier {
 private:

 public:
     CausalLiquefier();
     ~CausalLiquefier() {};

     void Init();
     void Exec();
     void Clear();

};

#endif  // CAUSALLIQUEFIER_H
