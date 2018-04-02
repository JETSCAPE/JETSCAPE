/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#ifndef COLORLESSHADRONIZATION_H
#define COLORLESSHADRONIZATION_H

#include "HadronizationModule.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

class ColorlessHadronization : public HadronizationModule<ColorlessHadronization>
{  
 public:

  ColorlessHadronization();
  virtual ~ColorlessHadronization();
  
  void Init();
  void DoHadronization(vector<vector<shared_ptr<Parton>>>& shower, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);

 protected:
  static Pythia8::Pythia pythia;  
  
};


#endif // COLORLESSHADRONIZATION_H

