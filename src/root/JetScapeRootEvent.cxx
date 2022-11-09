// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "JetScapeRootEvent.h"
#include "JetScapeLogger.h"

ClassImp(JetScapeRootEvent)

JetScapeRootEvent::JetScapeRootEvent()
{
  eventNumber=-99;
  VERBOSE(9);
}
  
JetScapeRootEvent::~JetScapeRootEvent()
{
  VERBOSE(9);
  //DEBUG;
}
