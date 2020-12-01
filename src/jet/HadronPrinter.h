/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2020
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


#ifndef HADRONPRINTER_H
#define HADRONPRINTER_H

#include "JetClass.h"
//#include <vector>
#include <string>
#include<fstream>
#include "JetScapeModuleBase.h"
//#include "PartonShower.h"
//#include "sigslot.h"

namespace Jetscape {

class HadronPrinter : public JetScapeModuleBase
{

 public:
 
  HadronPrinter();
  virtual ~HadronPrinter();
  
  virtual void Init();
  virtual void Exec() final;
  virtual void Clear();
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);

  sigslot::signal1<vector<shared_ptr<Hadron>>& > GetFinalHadronList;

  void SetFinalHadrons(vector<shared_ptr<Hadron>>& hadrons){
    finalHadrons = hadrons;
    PrintFinalHadron();
  }
  
  void PrintFinalHadron();
 
 private:

  vector<shared_ptr<Hadron>> finalHadrons;
	 std::ofstream fHadronOutfile;  ///< the output stream where events are saved to file

};

} // end namespace Jetscape


#endif
