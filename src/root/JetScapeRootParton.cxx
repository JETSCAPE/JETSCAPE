// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "JetScapeRootParton.h"
#include "JetScapeLogger.h"

ClassImp(JetScapeRootParton)

JetScapeRootParton::JetScapeRootParton() //: Parton()
{
  VERBOSE(9);
}
  
JetScapeRootParton::~JetScapeRootParton()
{
  VERBOSE(9);
  //DEBUG;
}
  
JetScapeRootParton::JetScapeRootParton(Parton& srp)
{
  label=srp.plabel();stat=srp.pstat();pid=srp.pid();
  pt=srp.pt();eta=srp.eta();phi=srp.phi();e=srp.e();
}
