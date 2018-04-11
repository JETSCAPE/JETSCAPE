/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * Modular, task-based framework
 * Intial Design: Joern Putschke, Kolja Kauder (Wayne State University)
 * For the full list of contributors see AUTHORS.

 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

/** TaskSupport instance class (meant as singelton)
 * Keeps track of every created task in a thread-safe manner
 * and provides resources that depend on such information
 * Initial reason was to provide tasks with random seeds
 * such that reproducible running is possible
 * It provides a random engine factory.
 * That way individual random engines can be used where necessary
 * but by default all can use the same 
 * (the state of mersenne twister can be  rather large)
 * 
 * Note that (apart from Microsoft VS2013), magic statics should ensure that
 * the Instance() method is automagically thread safe 
 * 
 * Note 2: make_unique doesn't work for some reason. "new" does the trick here though.
 */

#ifndef JETSCAPETASKSUPPORT_H
#define JETSCAPETASKSUPPORT_H

#include "InitialState.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "FluidDynamics.h"
#include "HardProcess.h"
#include "JetScapeWriter.h"

#include<iostream>
#include<atomic>
#include<memory>
#include<random>
#include<thread>

using std::atomic_int;

namespace Jetscape {

  class JetScapeTaskSupport //: public sigslot::has_slots<sigslot::multi_threaded_local>
  {
  
  public:

    static JetScapeTaskSupport* Instance();    
    // void CleanUp();

    /// Tasks should call this method at creation and
    /// remember the answer as their task id
    /// This could co a lot more, like keep a map of numbers to task.id (essentially the name of the task)
    /// But for now keep it simple
    int RegisterTask();

    /// Initialize random engine functionality from the XML file
    static void ReadSeedFromXML();

    /// Return a handle to a mersenne twister
    /// Usually this should be just one shared by everybody
    /// but if reproducible seeds are requested,
    /// every task gets their own
    shared_ptr<std::mt19937> get_mt19937_generator( int TaskId );
  
    // Getters
    static unsigned int GetRandomSeed() { return random_seed_;};

  protected:
    static bool one_generator_per_task_;
    
  private:

    JetScapeTaskSupport() : current_task_number_(0) {};

    static JetScapeTaskSupport* m_pInstance;
    
    atomic_int current_task_number_;
    static unsigned int random_seed_;
    static bool initialized_;
    
    static shared_ptr<std::mt19937> one_for_all_;
    
  };

} // end namespace Jetscape

#endif

