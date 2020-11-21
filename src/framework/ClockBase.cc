#include "ClockBase.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"

#include <iostream>

namespace Jetscape {

//bool ClockBase::use_clock = true;

ClockBase::ClockBase()
{
	id = "";
}

void ClockBase::Info()
{
    JSINFO<<GetId()<<" "<<GetTimeRefFrameId();
}

}