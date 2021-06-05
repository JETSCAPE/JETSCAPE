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

#ifndef TIMEMODULE_H
#define TIMEMODULE_H

#include "ModuleClock.h"
#include "MainClock.h"
//#include "JetScapeTask.h"

#include <string>
#include <memory>

using std::string;
using std::shared_ptr;

namespace Jetscape {

class TimeModule //: public JetScapeTask
{

public:

    TimeModule();

    TimeModule(double t1,double t2);

    virtual ~TimeModule() {};

    void ClockInfo();

    void AddModuleClock(shared_ptr<ModuleClock> m_mClock) {mClock = m_mClock;}
    shared_ptr<ModuleClock> GetModuleClock() const {return mClock;}

    void AddMainClock(shared_ptr<MainClock> m_mainClock);
    static shared_ptr<MainClock> GetMainClock() {return mainClock;}

    static bool ClockUsed() { return use_clock; }
    bool UseModuleClock() {if (mClock!=nullptr) return true; else return false;}

    //static bool use_clock; //better in time based module base ...

    double GetModuleCurrentTime();

    double GetModuleDeltaT();

    bool IsValidModuleTime() {if (GetModuleCurrentTime() >= t0 && GetModuleCurrentTime() < tn) return true; else return false;};

    void SetTimeRange(double t1, double t2) {t0 = t1; tn = t2;};

    const double GetTStart() const {return t0;};

    const double GetTEnd() const {return tn;};

private:

    shared_ptr<ModuleClock> mClock;
    static shared_ptr<MainClock> mainClock;

    double t0;// module start time; default is 0
    double tn;// module end time; default is 100

    static bool use_clock; //better in time based module base ...
};

} // end namespace Jetscape

#endif
