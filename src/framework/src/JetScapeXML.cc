// JetScape XML init reader class implementation (meant as singelton)

#include "JetScapeXML.h"
#include "JetScapeLogger.h"
#include <stdlib.h>

using namespace std;

JetScapeXML* JetScapeXML::m_pInstance = NULL;

JetScapeXML* JetScapeXML::Instance()
{
  if (!m_pInstance)
    {
      INFO<<"Created JetScapeXML Instance";
      m_pInstance = new JetScapeXML();
    }
  
  return m_pInstance;
}

void JetScapeXML::OpenXMLFile()
{
  if (!xml_file_open)
    {
      xml_doc.LoadFile((char*) GetXMLFileName().c_str());
      
      if (xml_doc.ErrorID()<1)
	{
	  // Also: __FILE__ in gcc & clang
	  //cout<< __PRETTY_FUNCTION__ <<":"<<__LINE__<< " Open XML file : "<< GetXMLFileName() << endl;
	  
	  INFO<<"Open XML file : "<< GetXMLFileName();
	  xml_root = (tinyxml2::XMLElement*) xml_doc.FirstChildElement("jetscape" );
	  
	  if (!xml_root)
	    {
	      WARN << "Not a valide JetScape XML file!";
	      exit(-1);
	    }
	}
      else
	{
	  WARN << "XML file not found/not properly opened! Error code : "<<xml_doc.ErrorID();
	  exit(-1);
	}
      xml_file_open=true;
    }
}

void JetScapeXML::OpenXMLFile(string m_name)
{
  SetXMLFileName(m_name);
  OpenXMLFile();
}
