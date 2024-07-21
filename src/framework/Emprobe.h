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
// JETSCAPE module for EM probe
// -----------------------------------------

#ifndef EMPROBE_H_
#define EMPROBE_H_

#include <vector>

#include "JetScapeModuleBase.h"
#include "JetClass.h"
#include "JetScapeWriter.h"
#include "FluidEvolutionHistory.h"

namespace Jetscape {

class Emprobe : public JetScapeModuleBase {
private:
  EvolutionHistory bulk_info;

public:
  Emprobe();
  ~Emprobe();

  virtual void Init();
  virtual void Exec();
  virtual void Clear();

  bool boost_invariance;
  bool check_boost_invariance();
};

} // end namespace Jetscape

#endif // SOFTPARTICLIZATION_H_
