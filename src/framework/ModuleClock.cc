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
	currentModuleTime = mainClock.lock()->GetCurrentTime() * 2.;
        moduleDeltaT = mainClock.lock()->GetDeltaT() * 2.;
}

}