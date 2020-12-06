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

//  Query History instance class (meant as singelton)

#ifndef QUERYHISTORY_H
#define QUERYHISTORY_H

#include "JetScapeModuleBase.h"
#include "cpp17/any.hpp"
#include "cpp17/variant.hpp"

#include <iostream>
#include <string>
#include <map>
#include <unordered_map>
//#include <cstddef>

//#include "sigslot.h"
//using namespace sigslot;

//JP: Make sure to not introduce and memory leak here using an instance.
//Should be fine (see SignalManager) by using weak pointers. But make sure !!!!!

//Maybe change, namespaces for any and varaint ...
using namespace linb;
using namespace mpark;

namespace Jetscape {

class QueryHistory
{
  public:

    static QueryHistory *Instance();
    
    void AddMainTask(std::shared_ptr<JetScapeTask> m_main_task) {main_task = m_main_task;}
    void UpdateTaskMap();
    void PrintTasks();
    void PrintTaskMap();

    std::unordered_multimap<std::string,std::weak_ptr<JetScapeTask> > GetTaskMap() {return taskMap;}

    //JP: same can be done with variant if all datatypes are know
    //and put into the varaint definition --> elevated to framework like data types
    //maybe not ideal, to be discussed ...
    any GetHistoryFromModule(string mName);

    //JP: maybe use as standard only to allow for multipe modules like in JetEnhergyLoss ...
    vector<any> GetHistoryFromModules(string mName);

  private:

    QueryHistory(){};
    QueryHistory(QueryHistory const &){};
    static QueryHistory *m_pInstance;

    std::unordered_multimap<std::string,std::weak_ptr<JetScapeTask> > taskMap;

    std::weak_ptr<JetScapeTask> main_task;

};

} // end namespace Jetscape

#endif