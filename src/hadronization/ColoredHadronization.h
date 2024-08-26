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

#ifndef COLOREDHADRONIZATION_H
#define COLOREDHADRONIZATION_H

#include "HadronizationModule.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

class ColoredHadronization : public HadronizationModule<ColoredHadronization> {
 public:
  ColoredHadronization();
  virtual ~ColoredHadronization();

  void Init();
  void DoHadronization(vector<vector<shared_ptr<Parton>>> &shower,
                       vector<shared_ptr<Hadron>> &hOut,
                       vector<shared_ptr<Parton>> &pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);

 private:
  double p_fake;

  // Allows the registration of the module so that it is available to be used by
  // the Jetscape framework.
  static RegisterJetScapeModule<ColoredHadronization> reg;

 protected:
  static Pythia8::Pythia pythia;
};

#endif  // COLOREDHADRONIZATION_H
