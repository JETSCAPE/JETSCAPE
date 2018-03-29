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

namespace Jetscape {

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

void JetScapeWriterAscii::WriteHeaderToFile()
{
  INFO<<"Run JetScapeWriterAscii: Write header of event # "<<GetCurrentEvent()<<" ...";
  Write(to_string(GetCurrentEvent()) + " Event");

  std::ostringstream oss;
  oss.str(""); oss << GetId() << "sigmaGen " << GetHeader().GetSigmaGen();
  WriteComment ( oss.str() );
  oss.str(""); oss << GetId() << "sigmaErr " << GetHeader().GetSigmaErr();
  WriteComment ( oss.str() );
  oss.str(""); oss << GetId() << "weight " << GetHeader().GetEventWeight();
  WriteComment ( oss.str() );
}
  
void JetScapeWriterAscii::WriteEvent()
{
  // INFO<<"Run JetScapeWriterAscii: Write event # "<<GetCurrentEvent()<<" ...";
  // do nothing, the modules handle this
}

void JetScapeWriterAscii::Write(weak_ptr<Parton> p)
{
  output_file<<*p.lock()<<endl;
}

void JetScapeWriterAscii::Write(weak_ptr<Vertex> v)
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
  // INFO<<"Run JetScapeWriterAscii: Write event # "<<GetCurrentEvent()<<" ...";
  
  // if (GetActive())
  //   WriteEvent();
}

void JetScapeWriterAscii::WriteInitFileXML()
{
  JSDEBUG<<"Write XML to output file. XML file = "<<JetScapeXML::Instance()->GetXMLFileName();
  tinyxml2::XMLPrinter printer;
  JetScapeXML::Instance()->GetXMLDocument().Print(&printer);
  WriteComment("Init XML file used : "+JetScapeXML::Instance()->GetXMLFileName());
  output_file<<printer.CStr();
}

void JetScapeWriterAscii::Write(weak_ptr<Hadron> h)
{
  output_file<<*h.lock()<<endl;
}

} // end namespace Jetscape
