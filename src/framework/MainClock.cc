#include "MainClock.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"

#include <iostream>

namespace Jetscape {

//bool ClockBase::use_clock = true;

MainClock::MainClock() : ClockBase()
{
    SetId("MainClock");
    deltaT = 0.1;
    startTime = 0.0;
    endTime = 20.0 + deltaT;
    currentTime = startTime;
}

MainClock::MainClock(string m_id, double m_st, double m_et, double m_dt) : ClockBase()
{
    SetId("MainClock");
    SetTimeRefFrameId(m_id); 
    deltaT = m_dt;
    startTime = m_st;
    endTime = m_et + deltaT;    
    currentTime = startTime;
}

MainClock& MainClock::operator++() {
    
    currentTime += deltaT;

    return *this;
  }

bool MainClock::Next()
{
    currentTime += deltaT;	
    if (currentTime<endTime)
	return true;
    else
	return false;
}

void MainClock::Info()
{
    ClockBase::Info();
    JSINFO<<" Start Time = "<<startTime<<" | End Time = "<<endTime;
    JSINFO<<" Curent Time = "<<currentTime;
}  

}