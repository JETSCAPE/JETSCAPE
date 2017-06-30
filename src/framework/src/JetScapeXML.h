// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// JetScape XML init reader class (meant as singelton)

#ifndef JETSCAPEXML_H
#define JETSCAPEXML_H

#include<iostream>
#include<string>

#include "tinyxml2.h"

using namespace std;

namespace Jetscape {

class JetScapeXML
{
  
 public:

  static JetScapeXML* Instance();

  tinyxml2::XMLElement* GetXMLRoot() {return xml_root;}
  tinyxml2::XMLDocument& GetXMLDocument() {return xml_doc;}

  void SetXMLFileName(string m_name) { xml_file_name = m_name; }
  string GetXMLFileName() {return xml_file_name;}
  bool IsOpen() {return xml_file_open;}
  
  void OpenXMLFile();
  void OpenXMLFile(string m_name);
  
 private:

  JetScapeXML() {xml_file_name="";xml_file_open=false;};
  JetScapeXML(JetScapeXML const&) {};
  static JetScapeXML* m_pInstance;

  tinyxml2::XMLElement *xml_root; //use unique pointer here instead of raw pointer (check with tinyxml interface)
  tinyxml2::XMLDocument xml_doc;

  string xml_file_name;
  bool xml_file_open;
  
};

} // end namespace Jetscape

#endif

