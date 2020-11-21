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
    endTime = 20.0;
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
    JSINFO<<"Start Time = "<<startTime<<" | End Time = "<<endTime;
    JSINFO<<"Curent Time = "<<currentTime;
}  

}