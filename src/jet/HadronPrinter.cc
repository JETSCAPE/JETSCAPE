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

#include "HadronPrinter.h"
namespace Jetscape {

HadronPrinter::HadronPrinter()
{
}

HadronPrinter::~HadronPrinter()
{
}

void HadronPrinter::Init()
{
  this->SetId("Printer");
  fHadronOutfile.open("finalStateHadrons.dat");
}

void HadronPrinter::Exec()
{
  VERBOSE(2) <<"Run HadronPrinter: ";
  GetFinalHadronList(finalHadrons);
  PrintFinalHadron();
}

void HadronPrinter::Clear()
  {this->finalHadrons.clear();}

void HadronPrinter::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);
  JetScapeTask::WriteTasks(w);
}

void HadronPrinter::PrintFinalHadron(){
  char buffer [33];

  JSINFO << "HadronPrinter : starting to print hadrons";
  int i=0;
    fHadronOutfile << "#  PID  pstat  E  Px  Py  Pz  eta  phi" << endl;
  for(auto it : finalHadrons){
    Hadron h = *it.get();
    fHadronOutfile << i <<"  "<< h.pid() <<"  "<< h.pstat() <<"  "<< h.e() <<"  "<< h.px() <<"  "<< h.py() <<"  "<< h.pz() <<"  "<< h.eta() <<"  "<< h.phi() << endl;
    ++i;
  }

}
} 


