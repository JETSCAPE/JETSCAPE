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

//REMARK: Old JetScape PSG w/o droplets etc ...
// pretty much copy of the DoShower in JetEnergyLoss ....

#ifndef PARTONSHOWERGENERATORDEFAULT_H
#define PARTONSHOWERGENERATORDEFAULT_H

#include "PartonShowerGenerator.h"

namespace Jetscape {

class JetEnergyLoss;

class PartonShowerGeneratorDefault : public PartonShowerGenerator
{
 public:

   PartonShowerGeneratorDefault() : PartonShowerGenerator()  {};
   virtual ~PartonShowerGeneratorDefault() {};

   virtual void DoShower(JetEnergyLoss &j);
};

} // end namespace Jetscape

#endif
