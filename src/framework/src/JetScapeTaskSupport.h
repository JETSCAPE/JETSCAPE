// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Kolja Kauder (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

/** TaskSupport instance class (meant as singelton)
 * Keeps track of every created task in a thread-safe manner
 * and provides resources that depend on such information
 * Initial reason was to provide tasks with random seeds
 * such that reproducible running is possible
 * Will eventually also have a random engine factory
 * That way individual random engines can be used where necessary
 * but by default all can use the same 
 * (the state of mersenne twister can be  rather large)
 * 
 * Note that (apart from Microsoft VS2013), magic statics should ensure that
 * the Instance() method is automagically thread safe 
 * 
 * Note 2: make_unique doesn't work for some reason. new does the trick here though.
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
#include<random>

 using namespace std;

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

    /// Return a handle to a mersenne twister
    /// Usually this should be just one shared by everybody
    /// (todo: what about threads?)
    /// but if reproducible seeds are requested,
    /// every task gets their own
  
  private:

  JetScapeTaskSupport(): CurrentTaskNumber(0) {};
    // JetScapeTaskSupport(JetScapeTaskSupport const&) {};
    static JetScapeTaskSupport* m_pInstance;

    atomic_int CurrentTaskNumber;
  
  };

} // end namespace Jetscape

#endif

