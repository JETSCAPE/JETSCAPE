// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// jetscape writer ascii class

#include "JetScapeWriterAscii.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"

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
  DEBUG<< GetCurrentEvent() << " Event";
  Write(to_string(GetCurrentEvent()) + " Event");
  
}

void JetScapeWriterAscii::Write(weak_ptr<Parton> p)
{
  output_file<<*p.lock()<<endl;
}

void JetScapeWriterAscii::Write(weak_ptr<VertexBase> v)
{
  output_file<<*v.lock()<<endl;
}

void JetScapeWriterAscii::Init()
{
   if (GetActive())
     {
       INFO<<"JetScape Ascii Writer initialized with output file = "<<GetOutputFileName();
       output_file.open(GetOutputFileName().c_str());
       
       //Write Init Informations, like XML and ... to file ...
       //WriteInitFileXML();
     }
}

void JetScapeWriterAscii::Exec()
{
  if (GetActive())
    WriteEvent();
}

void JetScapeWriterAscii::WriteInitFileXML()
{
  DEBUG<<"Write XML to output file. XML file = "<<JetScapeXML::Instance()->GetXMLFileName();
  tinyxml2::XMLPrinter printer;
  JetScapeXML::Instance()->GetXMLDocument().Print(&printer);
  //cout<<printer.CStr()<<endl;
  WriteComment("Init XML file used : "+JetScapeXML::Instance()->GetXMLFileName());
  output_file<<printer.CStr();
}
