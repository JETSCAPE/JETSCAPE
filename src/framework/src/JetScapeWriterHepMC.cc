// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "JetScapeWriterHepMC.h"
#include "JetScapeLogger.h"


JetScapeWriterHepMC::~JetScapeWriterHepMC()
{
  if (GetActive())
      Close();
}

void JetScapeWriterHepMC::WriteEvent()
{
  INFO<< GetCurrentEvent() << " in HepMC ... ";
  
  GenEvent evt(Units::GEV,Units::MM);
  Print::content(evt);
  write_event(evt);
}

void JetScapeWriterHepMC::Init()
{
   if (GetActive())
     {
       INFO<<"JetScape HepMC Writer initialized with output file = "<<GetOutputFileName();
     }
}

void JetScapeWriterHepMC::Exec()
{
  if (GetActive())
    WriteEvent();
}
