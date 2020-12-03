#include "TimeModule.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"

#include <iostream>

namespace Jetscape {

bool TimeModule::use_clock = false;
shared_ptr<MainClock> TimeModule::mainClock = nullptr;

TimeModule::TimeModule() 
{
    mClock = nullptr;    
}

void TimeModule::ClockInfo()
{
    JSINFO<<"TimeModule::Info()";
    JSINFO<<"Main clocked used = "<<ClockUsed();
    JSINFO<<"Module clock = "<<UseModuleClock();	
    if (UseModuleClock())
    	mClock->Info();
    else if (ClockUsed())
    	mainClock->Info();
    JSINFO<<"Current Module Time = "<<GetModuleCurrentTime();
}

double TimeModule::GetModuleDeltaT()
{
    double dT=-99.0;

    if (UseModuleClock())
    {
        mClock->Transform(mainClock);
        dT = mClock->GetDeltaT();
    }
    else if (ClockUsed())
    {
        dT = mainClock->GetDeltaT();
    }
    else
        {JSWARN<<"No clocks found ..."; exit(-1);}

    return dT;
}

double TimeModule::GetModuleCurrentTime()
{
    double cTime=-99.0;

    if (UseModuleClock())
    {
    	mClock->Transform(mainClock);
        cTime = mClock->GetCurrentTime();
    }
    else if (ClockUsed())
    {
        cTime = mainClock->GetCurrentTime();
    }
    else
    	{JSWARN<<"No clocks found ..."; exit(-1);}

    return cTime;
}

}