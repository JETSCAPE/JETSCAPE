/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/amajumder/JETSCAPE-COMP/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#ifndef PARTONSHOWERGENERATOR_H
#define PARTONSHOWERGENERATOR_H

namespace Jetscape {

class JetEnergyLoss;

class PartonShowerGenerator
{
 public:

  PartonShowerGenerator() {};
  virtual ~PartonShowerGenerator() {};

  virtual void DoShower(JetEnergyLoss &j);


};

} // end namespace Jetscape

#endif
