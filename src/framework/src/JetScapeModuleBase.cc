// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "JetScapeModuleBase.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"

#include<iostream>

namespace Jetscape {

int JetScapeModuleBase::current_event=0;

JetScapeModuleBase::JetScapeModuleBase()
{
  xml_file_name = "";
}

JetScapeModuleBase::JetScapeModuleBase(string m_name)
{
  xml_file_name = m_name;
}

JetScapeModuleBase::~JetScapeModuleBase()
{
  disconnect_all();
}

void JetScapeModuleBase::Init()
{
  if (!JetScapeXML::Instance()->GetXMLRoot())
    {
      WARN << "Not a valide JetScape XML file or no XML file loaded!";
      exit(-1);
    }
}

} // end namespace Jetscape
