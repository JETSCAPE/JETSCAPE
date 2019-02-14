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

#include "JetScape.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include "JetEnergyLossManager.h"
#include "FluidDynamics.h"
#include "JetScapeBanner.h"
#include "InitialState.h"
#include "PreequilibriumDynamics.h"

#include<iostream>

using namespace std;

namespace Jetscape {

  /** Default constructor to create the main task of the JetScape framework. It sets the total number of events to 1.
   * By default, hydro events are used only once
   */
JetScape::JetScape() 
  : n_events(1)
  , reuse_hydro_ (false)
  , n_reuse_hydro_ (0)
{
  VERBOSE(8);
}

JetScape::~JetScape()
{
  VERBOSE(8);
  //JetScapeSignalManager::Instance()->Clear();
  //not needed, use weak_ptr in JetScapeSignalManager class (=not owning)
}

void JetScape::Show()
{
  ShowJetscapeBanner();
}

void JetScape::Init()
{
  Show();
  
  JSINFO<<BOLDRED<<"Intialize JetScape ...";
  
  JetScapeXML::Instance()->OpenXMLFile(GetXMLFileName());
  
  string log_debug = JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement( "debug" )->GetText();
  string log_remark = JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement( "remark" )->GetText();
  int m_vlevel=0;
  
  JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("vlevel")->QueryIntText(&m_vlevel);
  
  // Some general JetScape interface settings can also be read in from XML file ...
  // Order of task could also be verifed via XML or via sort of TaskList with order parameter
  if ((int) log_debug.find("off")>=0)
    JetScapeLogger::Instance()->SetDebug(false);
  
  if ((int) log_remark.find("on")>=0)
    JetScapeLogger::Instance()->SetRemark(true);

  JSDEBUG<<"JetScape Debug from XML = "<< log_debug;
  JSDEBUG<<"JetScape Remark from XML = "<< log_remark;

   if ((int) m_vlevel>0)
     {
       JetScapeLogger::Instance()->SetVerboseLevel(m_vlevel);
       JSDEBUG<<"JetScape Verbose Level from XML = "<<m_vlevel;
     }

   SetPointers();
   
   // Set up helper. Mostly used for random numbers
   // Needs the XML reader singleton set up
   JSDEBUG<<"Seeding JetScapeTaskSupport from XML";
   JetScapeTaskSupport::ReadSeedFromXML( );

  // Has to be called explicitly since not really fully recursively (if ever needed)
   // So --> JetScape is "Task Manager" of all modules ...
   
   JSINFO<<"Found "<<GetNumberOfTasks()<<" Modules Initialize them ... ";
   JetScapeTask::InitTasks();
}

// kind of cluncky, maybe a better way ... ?
// Handle signal/slots in JetScape hence avoid passing pointers to sub tasks ...
void JetScape::SetPointers() {
    // to get hydro pointer for signals, use signal?
    JSINFO << "Set Hydro,JetEnergylossManager and IS Pointers for SignalManager "
         << "to create Signal/Slots";
    
    bool hydro_pointer_is_set = false;
    for (auto it : GetTaskList()) {
        if (dynamic_pointer_cast<InitialState>(it)) {
            JetScapeSignalManager::Instance()->SetInitialStatePointer(
                                    dynamic_pointer_cast<InitialState>(it));
        } else if (dynamic_pointer_cast<PreequilibriumDynamics>(it)) {
            JetScapeSignalManager::Instance()->SetPreEquilibriumPointer(
                        dynamic_pointer_cast<PreequilibriumDynamics>(it));
        } else if (dynamic_pointer_cast<FluidDynamics>(it)
                   && !hydro_pointer_is_set) {
            JetScapeSignalManager::Instance()->SetHydroPointer(
                        dynamic_pointer_cast<FluidDynamics>(it));
            hydro_pointer_is_set = true;
        } else if (dynamic_pointer_cast<JetEnergyLossManager>(it)) {
            JetScapeSignalManager::Instance()->SetJetEnergyLossManagerPointer(
                        dynamic_pointer_cast<JetEnergyLossManager>(it));
        } else if (dynamic_pointer_cast<HardProcess>(it)) {
            JetScapeSignalManager::Instance()->SetHardProcessPointer(
                        dynamic_pointer_cast<HardProcess>(it));
        } else if (dynamic_pointer_cast<JetScapeWriter>(it)
                   && it->GetActive()) {
            JetScapeSignalManager::Instance()->SetWriterPointer(
                        dynamic_pointer_cast<JetScapeWriter>(it)); 
        } else if (dynamic_pointer_cast<PartonPrinter>(it)) {
            JetScapeSignalManager::Instance()->SetPartonPrinterPointer(
                        dynamic_pointer_cast<PartonPrinter>(it));
        } else if (dynamic_pointer_cast<SoftParticlization>(it)) {
            JetScapeSignalManager::Instance()->SetSoftParticlizationPointer(
                        dynamic_pointer_cast<SoftParticlization>(it));
        }
    }
}

void JetScape::Exec()
{
  JSINFO<<BOLDRED<<"Run JetScape ...";
  JSINFO<<BOLDRED<<"Number of Events = "<<GetNumberOfEvents();
  
  // JetScapeTask::ExecuteTasks(); Has to be called explicitly since not really fully recursively (if ever needed)
  // --> JetScape is "Task Manager" of all modules ...

  // Simple way of passing the writer module pointer
  vector<weak_ptr<JetScapeWriter>> vWriter;
  
  for (auto it : GetTaskList()) {
    if (dynamic_pointer_cast<JetScapeWriter>(it)){  
      if (it->GetActive())
	vWriter.push_back(dynamic_pointer_cast<JetScapeWriter>(it));
    }
  } 
  
  for (int i=0;i<GetNumberOfEvents();i++)
    {
      JSINFO<<BOLDRED<<"Run Event # = "<<i;
      JSDEBUG<<"Found "<<GetNumberOfTasks()<<" Modules Execute them ... ";

      // First run all tasks
      JetScapeTask::ExecuteTasks();

      // Then hand around the collection of writers and ask
      // modules to write what they like
      // Sequence of events:
      // -- writer->Exec is called and redirects to WriteEvent, which starts a new event line
      // -- any remaining exec's finish
      // -- all modules write their headers
      // -- Now all header info is known to the writers, so write out the header
      // -- all other Write()'s are being called
      // the result still confuses me. It's in the best possible order but it shouldn't be.
      
      // collect module header data
      for (auto w : vWriter) {
	auto f = w.lock();
	if ( f ) JetScapeTask::CollectHeaders(w);
      }
      // official header
      for (auto w : vWriter) {
	auto f = w.lock();
	if ( f ) f->WriteHeaderToFile();
      }
      // event data
      for (auto w : vWriter) {
	auto f = w.lock();
	if ( f ) JetScapeTask::WriteTasks(w);
      }

      // Finalize
      for (auto w : vWriter) {
	auto f = w.lock();
	if ( f ) f->WriteEvent();
      }

      // For reusal, deactivate task after it has finished but before it gets cleaned up.
      if (reuse_hydro_) {
	if (n_reuse_hydro_ <= 0) {
	  JSWARN << " reuse_hydro is set, but n_reuse_hydro=" << n_reuse_hydro_;
	  throw std::runtime_error ("Incompatible reusal settings.");
	}
	for (auto it : GetTaskList()) {
	  if (!dynamic_pointer_cast<FluidDynamics>(it)
	      && !dynamic_pointer_cast<InitialState>(it)) {
	    continue;
	  }
	  if (i%n_reuse_hydro_ == n_reuse_hydro_ - 1) {
	    JSDEBUG << " i was " << i << " i%n_reuse_hydro_ = " << i%n_reuse_hydro_ << " --> ACTIVATING";
	    it->SetActive(true);
	  } else {
	    JSDEBUG << " i was " << i << " i%n_reuse_hydro_ = " << i%n_reuse_hydro_ << " --> DE-ACTIVATING";
	    it->SetActive(false);
	  }
	}
      }
      // Now clean up, only affects active taskjs
      JetScapeTask::ClearTasks();

      IncrementCurrentEvent();
    }
}

void JetScape::Finish()
{
  JSINFO<<BOLDBLACK<<"JetScape finished after "<<GetNumberOfEvents()<<" events!";
  JSDEBUG<<"More infos wrap up/saving to file/closing file ...";

  // same as in Init() and Exec() ...
  JetScapeTask::FinishTasks(); //dummy so far ...
}

} // end namespace Jetscape
