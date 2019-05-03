/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#ifndef JETSCAPEXML_H
#define JETSCAPEXML_H

#include<iostream>
#include<string>
#include<stdexcept>
#include<initializer_list>

#include "tinyxml2.h"

/**
 * @class JetScapeXML
 * @brief JetScape XML init reader class (meant as singleton)
 *
 * This class contains the machinery to load two XML configuration files: a Master file, and a User file.
 *
 */

using std::string;
using std::runtime_error;

namespace Jetscape {

class JetScapeXML
{
  
 public:

  static JetScapeXML* Instance();
  
  // Master file
  
  tinyxml2::XMLElement* GetXMLRootMaster() {return xml_root_master;}
  tinyxml2::XMLDocument& GetXMLDocumentMaster() {return xml_doc_master;}
  tinyxml2::XMLElement* GetXMLElementMaster(std::initializer_list<const char*> &path);

  void SetXMLMasterFileName(string m_name) { xml_master_file_name = m_name; }
  std::string GetXMLMasterFileName() {return xml_master_file_name;}
  bool IsMasterFileOpen() {return xml_master_file_open;}
  
  void OpenXMLMasterFile();
  void OpenXMLMasterFile(string m_name);
  
  // User file

  tinyxml2::XMLElement* GetXMLRootUser() {return xml_root_user;}
  tinyxml2::XMLDocument& GetXMLDocumentUser() {return xml_doc_user;}
  tinyxml2::XMLElement* GetXMLElementUser(std::initializer_list<const char*> &path);

  void SetXMLUserFileName(string m_name) { xml_user_file_name = m_name; }
  std::string GetXMLUserFileName() {return xml_user_file_name;}
  bool IsUserFileOpen() {return xml_user_file_open;}
  
  void OpenXMLUserFile();
  void OpenXMLUserFile(string m_name);
  
  // Helper functions for XML parsing/
  // Look first in user XML file for a parameter, and if not found look in the master XML file.
  tinyxml2::XMLElement* GetElement(std::initializer_list<const char*> path, bool isRequired = true);
  std::string GetElementText(std::initializer_list<const char*> path, bool isRequired = true);
  int GetElementInt(std::initializer_list<const char*> path, bool isRequired = true);
  double GetElementDouble(std::initializer_list<const char*> path, bool isRequired = true);
  
 private:
  
  JetScapeXML() {xml_master_file_name="";xml_master_file_open=false; xml_user_file_name="";xml_user_file_open=false;};
  JetScapeXML(JetScapeXML const&) {};
  static JetScapeXML* m_pInstance;

  // Master file
  
  tinyxml2::XMLElement *xml_root_master; //use unique pointer here instead of raw pointer (check with tinyxml interface)
  tinyxml2::XMLDocument xml_doc_master;
  
  std::string xml_master_file_name;
  bool xml_master_file_open;
  
  // User file

  tinyxml2::XMLElement *xml_root_user; //use unique pointer here instead of raw pointer (check with tinyxml interface)
  tinyxml2::XMLDocument xml_doc_user;

  std::string xml_user_file_name;
  bool xml_user_file_open;
  
};
  
// Print the XML element path name
std::ostream& operator << (std::ostream& os, std::initializer_list<const char*> path);

} // end namespace Jetscape

#endif

