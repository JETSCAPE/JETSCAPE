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

//REMARK JP: Current module transform for testing, just *2 !!!

#ifndef MODULECLOCK_H
#define MODULECLOCK_H

#include "ClockBase.h"
#include "MainClock.h"
#include <string>
#include <memory>

using std::string;

namespace Jetscape {

class ModuleClock 
	: public ClockBase {

public:

	ModuleClock();
	virtual ~ModuleClock() {};

	//virtual void Transform(string mainClockRef, double mainClockCurrentTime);
	virtual void Transform(std::weak_ptr<MainClock> mainClock);
	virtual void Info();

	inline double GetCurrentTime() {return currentModuleTime;}
        inline double GetDeltaT() {return moduleDeltaT;}
private:

	double currentModuleTime;
        double moduleDeltaT;

};

} // end namespace Jetscape

#endif