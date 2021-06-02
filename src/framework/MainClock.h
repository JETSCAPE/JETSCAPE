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

// Remark JP: Just basic functionalities, has to be extended and potentially made a bit smarter ...
// Currently assuming uniform ... but should be able to be overwritten via Next and ++ operator
// Include reading parameters from XML ... (to be done ...)

#ifndef MAINCLOCK_H
#define MAINCLOCK_H

#include "ClockBase.h"
#include <string>
#include <memory>

using std::string;

namespace Jetscape {

class MainClock
    : public ClockBase {

public:

    MainClock();
    MainClock(string m_id, double m_st, double m_et, double m_dt);

    virtual ~MainClock() {};

    void Reset() {currentTime = startTime;}
    void ResetToTime(double m_resetTime) {currentTime = m_resetTime;}

    virtual MainClock& operator++(); //check usage ... maybe better to always use Next() function ...
    virtual bool Next(); //check if max time check usage sufficent ...
    virtual bool Tick() {return Next();} 
    void Info();

    void SetStartTime(double m_StartTime) {startTime = m_StartTime;}
    void SetEndTime(double m_EndTime) {endTime = m_EndTime;}
    void SetDeltaT(double m_deltaT) {deltaT = m_deltaT;}

    inline double GetStartTime() {return startTime;}
    inline double GetEndTime() {return endTime;}
    inline double GetDeltaT() {return deltaT;}

    inline double GetCurrentTime() {return currentTime;}

private:

    double startTime;
    double endTime;
    double deltaT;
    double currentTime;

};

} // end namespace Jetscape

#endif