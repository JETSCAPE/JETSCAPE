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

#ifndef PARTONSHOWERGENERATOR_H
#define PARTONSHOWERGENERATOR_H

namespace Jetscape {

class JetEnergyLoss;

class PartonShowerGenerator {

public:

  PartonShowerGenerator(){};
  virtual ~PartonShowerGenerator(){};

  virtual void DoShower(JetEnergyLoss &j) {};

  virtual void DoCalculateTime(JetEnergyLoss &j) {};
  virtual void DoExecTime(JetEnergyLoss &j) {};
  virtual void DoInitPerEvent(JetEnergyLoss &j) {};
  virtual void DoFinishPerEvent(JetEnergyLoss &j) {};
  
};

} // end namespace Jetscape

#endif
