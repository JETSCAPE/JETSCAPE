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
// This is a wrapper for 
// -----------------------------------------

#ifndef EMPROBEWRAPPER_H
#define EMPROBEWRAPPER_H

#include <memory>

#include "Emprobe.h"
#include "JS_dilepton.h"

using namespace Jetscape;

class DileptonWrapper : public Emprobe {
private:
  tinyxml2::XMLElement *emprobe_xml_;

  int statusCode_;
  std::unique_ptr<Dilepton::JS_dilepton> Dilepton_ptr_;

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<DileptonWrapper> reg;
  vector<float> bulk_info_array;
  std::shared_ptr<std::vector<float>> dilepton_spec;  

public:
  DileptonWrapper();
  ~DileptonWrapper();

  void InitTask();
  void Exec();
  void Clear();
  void WriteTask(weak_ptr<JetScapeWriter> w);
  void getBulkInforfromJetScape();
};

#endif // 
