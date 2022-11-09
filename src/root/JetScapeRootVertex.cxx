// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "JetScapeRootVertex.h"
#include "JetScapeLogger.h"

ClassImp(JetScapeRootVertex)

JetScapeRootVertex::JetScapeRootVertex()
{
  VERBOSE(9);
}
  
JetScapeRootVertex::~JetScapeRootVertex()
{
  VERBOSE(9);
}
  
JetScapeRootVertex::JetScapeRootVertex(Vertex& srp)
{
  x=srp.x_in().x();y=srp.x_in().y();z=srp.x_in().z();t=srp.x_in().t();
}
