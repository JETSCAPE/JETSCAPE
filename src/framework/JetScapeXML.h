/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/amajumder/JETSCAPE-COMP/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

// JetScape XML init reader class (meant as singelton)

#ifndef JETSCAPEXML_H
#define JETSCAPEXML_H

#include<iostream>
#include<string>
#include<stdexcept>

#include "tinyxml2.h"

using std::string;
using std::runtime_error;

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

