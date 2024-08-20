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

#ifndef JETSCAPE_H
#define JETSCAPE_H

#if defined(__linux__)
#include <sys/sysinfo.h>
#include <sys/utsname.h>
#include <unistd.h>
#endif
#include "JetScapeLogger.h"
#include "JetScapeTaskSupport.h"
#include "JetScapeModuleBase.h"
#include "CausalLiquefier.h"

namespace Jetscape {

class JetScape : public JetScapeModuleBase {

public:
  /** Default constructor to create the main task of the JetScape framework.
  */
  JetScape();

  /** This is a destructor for a JetScape.
   */
  virtual ~JetScape();

  /** This function initializes the main task of the JetScape framework. It calls JetScapeTask::InitTaks() function to initialize the modules/tasks of a JetScapeTask.
  */
  void Init();

  /** This function execute the modules/tasks of a JetScapeTask for all the events. It also calls "GetPartons()" function to print parton shower, and  "WriteTasks()" function to store the data in the XML file.  
  */
  void Exec();

  void Finish();

  /** This function sets the total number of events to "m_n_events".
   */
  void SetNumberOfEvents(int m_n_events) { n_events = m_n_events; }

  /** This function returns the total number of events.
   */
  int GetNumberOfEvents() { return n_events; }

  /** Controls whether to reuse a hydro event (for speedup).
      The number of times is controled by SetNReuseHydro
   */
  inline void SetReuseHydro(const bool reuse_hydro) {
    reuse_hydro_ = reuse_hydro;
  }
  /** Returns whether hydro events are reused.
   */
  inline bool GetReuseHydro() const { return reuse_hydro_; }

  /** Controls number of times a hydro event gets reused.
      Reusal has to be explicitly turned on by SetReuseHydro.
      Turn it on first to avoid a warning.
   */
  inline void SetNReuseHydro(const unsigned int n_reuse_hydro) {
    if (!GetReuseHydro()) {
      JSWARN << "Number of hydro reusals set, but reusal not turned on.";
      JSWARN << "Try jetscape->SetReuseHydro (true);";
    }
    n_reuse_hydro_ = n_reuse_hydro;
  }
  inline unsigned int GetNReuseHydro() const { return n_reuse_hydro_; }

protected:
#if defined(__linux__)
  void PrintHardwareInfo();
#endif
  void CompareElementsFromXML();
  void recurseToBuild(std::vector<std::string> &elems, tinyxml2::XMLElement *mElement);
  void recurseToSearch(std::vector<std::string> &elems, tinyxml2::XMLElement *uElement);
  void ReadGeneralParametersFromXML();
  void DetermineTaskListFromXML();
  void DetermineWritersFromXML();
  void CheckForWriterFromXML(const char *writerName,
                             std::string outputFilename);
  void SetModuleId(tinyxml2::XMLElement *moduleElement,
                   shared_ptr<JetScapeModuleBase> module);

  void SetPointers();

  void Show();
  int n_events;
  int n_events_printout;

  bool reuse_hydro_;
  unsigned int n_reuse_hydro_;

  std::shared_ptr<CausalLiquefier> liquefier;

  bool
      fEnableAutomaticTaskListDetermination; // Option to automatically determine the task list from the XML file,
      // rather than manually calling JetScapeTask::Add() in the run macro.
};

} // end namespace Jetscape

#endif
