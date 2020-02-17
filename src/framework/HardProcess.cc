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

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include <string>

#include<iostream>

using namespace std;

#define MAGENTA "\033[35m"

namespace Jetscape {

HardProcess::HardProcess()
{
  VERBOSE(8);
  SetId("HardProcess");
}

HardProcess::~HardProcess()
{
  VERBOSE(8);
  hp_list.clear();
  hd_list.clear();
  disconnect_all();
}

void HardProcess::Init()
{
  JetScapeModuleBase::Init();

  JSINFO<<"Intialize HardProcess : "<<GetId()<< " ...";
  
  VERBOSE(8);

  ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
  if (!ini) {
    
    // If not vacuum case, give warning to add initial state module
    bool in_vac = GetXMLElementInt({"Eloss", "Matter", "in_vac"});
    if (!in_vac) {
      JSWARN << "No initial state module! Please check whether you intend to add an initial state module.";
    }
  }
  
  InitTask();

  JetScapeTask::InitTasks();
}

void HardProcess::Exec()
{
  JSINFO<<"Run Hard Process : "<<GetId()<< " ...";
  VERBOSE(8)<<"Current Event #"<<GetCurrentEvent();
  
  JetScapeTask::ExecuteTasks();
}

void HardProcess::Clear()
{
  JSDEBUG<<"Clear Hard Process : "<<GetId()<< " ...";

  hp_list.clear();
  hd_list.clear();
  VERBOSE(8)<<hp_list.size();
}

void HardProcess::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);

  auto f = w.lock();
  if ( f ){
    VERBOSE(8)<<f->GetOutputFileName();
    
    // Weight, xsec, etc
 
    // // Can explicitly write our own header information, though the writer should handle this.
    // std::ostringstream oss;
    // oss.str(""); oss << GetId() << " sigmaGen  = " << GetSigmaGen();  
    // f->WriteComment ( oss.str() );
    // oss.str(""); oss << GetId() << " sigmaErr  = " << GetSigmaErr();
    // f->WriteComment ( oss.str() );
    // oss.str(""); oss << GetId() << " weight  = " << GetEventWeight();
    // f->WriteComment ( oss.str() );
    
    // Hard partons
    f->WriteComment("HardProcess Parton List: "+GetId());  
    for ( auto hp : hp_list )    f->Write( hp );
  }
}

void HardProcess::CollectHeader( weak_ptr<JetScapeWriter> w ){
  auto f = w.lock();
  if ( f ){
    auto& header = f->GetHeader();
    header.SetSigmaGen( GetSigmaGen() );
    header.SetSigmaErr( GetSigmaErr() );
    header.SetEventWeight( GetEventWeight() );
  }
}

    
} // end namespace Jetscape
