// jetscape writer ascii + gzip class implementation

#include "JetScapeWriterAscii.h"
#include "JetScapeLogger.h"

JetScapeWriterAscii::JetScapeWriterAscii(string m_file_name_out)
{
  SetOutputFileName(m_file_name_out);
}

JetScapeWriterAscii::~JetScapeWriterAscii()
{
  VERBOSE(8);
  if (GetActive())
      Close();
}

void JetScapeWriterAscii::WriteEvent()
{
  DEBUG<< GetCurrentEvent() << " in Ascii ... ";
  Write(to_string(GetCurrentEvent()) + " in Ascii ... ");
  
}

void JetScapeWriterAscii::Write(shared_ptr<Parton> p)
{
  output_file<<*p<<endl;
}

void JetScapeWriterAscii::Init()
{
   if (GetActive())
     {
       INFO<<"JetScape Ascii Writer initialized with output file = "<<GetOutputFileName();
       output_file.open(GetOutputFileName().c_str());
     }
}

void JetScapeWriterAscii::Exec()
{
  if (GetActive())
    WriteEvent();
}
