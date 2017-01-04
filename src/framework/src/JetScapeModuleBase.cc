// JetScapeModuleBase class implementation
#include "JetScapeModuleBase.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"

#include<iostream>

int JetScapeModuleBase::current_event=0;

JetScapeModuleBase::JetScapeModuleBase()
{
  //Simple Debug replace --> logger
  //cout<<"JetScapeModuleBase : Default Constructor called."<<endl;
  xml_file_name = "";
}

JetScapeModuleBase::JetScapeModuleBase(string m_name)
{
  xml_file_name = m_name;
}

JetScapeModuleBase::~JetScapeModuleBase()
{
  //Simple Debug replace --> logger
  //cout<<"JetScapeModuleBase : Default Destructor called."<<endl;
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
