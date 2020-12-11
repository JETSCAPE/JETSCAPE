#include "ModuleClock.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"

#include <iostream>

namespace Jetscape {

//bool ClockBase::use_clock = true;

ModuleClock::ModuleClock() : ClockBase()
{
	SetId("ModuleClock");
	currentModuleTime = -99.;
        moduleDeltaT = -99.;
}

void ModuleClock::Info()
{
	ClockBase::Info();
	JSINFO<<" Curent Module Time = "<<currentModuleTime;
}  

/*
void ModuleClock::Transform(string mainClockRef, double mainClockCurrentTime)
{
	currentModuleTime = mainClockCurrentTime * 2.;
}
*/

void ModuleClock::Transform(std::weak_ptr<MainClock> mainClock)
{
	auto p = mainClock.lock();
	if (p)
	{
		currentModuleTime = p->GetCurrentTime() * 2.;
		moduleDeltaT = p->GetDeltaT() * 2.;
	}
	else
	{JSWARN << "Trying to transform module clock with no main clock ... exiting ..."; exit(-1);}

}
    
}
