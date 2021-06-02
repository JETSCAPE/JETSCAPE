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
TimeModule::TimeModule(double t1, double t2)
{
  mClock = nullptr;
  t0 = t1;
  tn = t2;
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
    
void TimeModule::AddMainClock(shared_ptr<MainClock> m_mainClock)
{
    if(mainClock != nullptr)
	JSWARN<<"Trying to add more than one main clock - will keep first one only ...";
    else
	mainClock = m_mainClock;
    use_clock = true;
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
