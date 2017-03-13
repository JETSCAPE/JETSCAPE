// jetscape writer ascii + gzip class implementation

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
