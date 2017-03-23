// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// jetscape writer ascii + gzip class implementation

#include "JetScapeWriterAsciiGZ.h"
#include "JetScapeLogger.h"

JetScapeWriterAsciiGZ::JetScapeWriterAsciiGZ(string m_file_name_out)
{
  SetOutputFileName(m_file_name_out);
}

JetScapeWriterAsciiGZ::~JetScapeWriterAsciiGZ()
{
  if (GetActive())
      Close();
}

void JetScapeWriterAsciiGZ::WriteEvent()
{
  DEBUG<< GetCurrentEvent() << " in Ascii gzip ... ";
  output_file<< GetCurrentEvent() << " in Ascii gzip ... \n";
}

void JetScapeWriterAsciiGZ::Write(shared_ptr<Parton> p)
{
  output_file<<*p<<endl;
}

void JetScapeWriterAsciiGZ::Init()
{
   if (GetActive())
     {
       INFO<<"JetScape (gzip) Writer initialized with output file = "<<GetOutputFileName();
       output_file.open(GetOutputFileName().c_str());
     }
}

void JetScapeWriterAsciiGZ::Exec()
{
  if (GetActive())
    WriteEvent();
}
