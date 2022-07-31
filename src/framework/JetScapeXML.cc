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

#include "JetScapeXML.h"
#include "JetScapeLogger.h"
#include <stdlib.h>

using namespace std;

namespace Jetscape {

JetScapeXML *JetScapeXML::m_pInstance = NULL;

JetScapeXML *JetScapeXML::Instance() {
  if (!m_pInstance) {
    JSINFO << "Created JetScapeXML Instance";
    m_pInstance = new JetScapeXML();
  }

  return m_pInstance;
}

//________________________________________________________________
void JetScapeXML::OpenXMLMasterFile() {
  JSWARN << "Deprecated function OpenXMLMasterFile(): Call OpenXMLMainFile() instead!";
  OpenXMLMainFile();
}

void JetScapeXML::OpenXMLMainFile() {

  if (!IsMainFileOpen()) {

    xml_doc_main.LoadFile((char *)GetXMLMainFileName().c_str());
    VERBOSE(2) << "Trying XML Main file : " << GetXMLMainFileName();

    if (xml_doc_main.ErrorID() < 1) {
      JSINFO << "Open XML Main file : " << GetXMLMainFileName();
      xml_root_main =
          (tinyxml2::XMLElement *)xml_doc_main.FirstChildElement("jetscape");

      if (!xml_root_main) {
        JSWARN << "Not a valid JetScape XML Main file!";
        exit(-1);
      }
    } else { //Check for an old default Master file
	  auto errCode = xml_doc_main.ErrorID(); //Save the original error code
      SetXMLMainFileName("../config/jetscape_master.xml"); //Try old default Master XML file
      xml_doc_main.LoadFile((char *)GetXMLMainFileName().c_str());
      VERBOSE(2) << "Looking for Old XML Master file : " << GetXMLMainFileName();

      if (xml_doc_main.ErrorID() < 1) {
        JSWARN << "Using Deprecated XML Master file : " << GetXMLMainFileName();
        xml_root_main =
            (tinyxml2::XMLElement *)xml_doc_main.FirstChildElement("jetscape");

        if (!xml_root_main) {
          JSWARN << "Not a valid JetScape XML Main file!";
          exit(-1);
        }
      } else {
          JSWARN << "XML Main file not found/not properly opened! Error code : "
                 << errCode;
          exit(-1);
	      }
      }

    xml_main_file_open = true;
  }
}

//________________________________________________________________
void JetScapeXML::OpenXMLMasterFile(string m_name) {
  JSWARN << "Deprecated function OpenXMLMasterFile(): Call OpenXMLMainFile() instead!";
  OpenXMLMainFile(m_name);
}

//________________________________________________________________
void JetScapeXML::OpenXMLMainFile(string m_name) {
  SetXMLMainFileName(m_name);
  OpenXMLMainFile();
}

//________________________________________________________________
void JetScapeXML::OpenXMLUserFile() {

  if (!xml_user_file_open) {

    xml_doc_user.LoadFile((char *)GetXMLUserFileName().c_str());
    VERBOSE(2) << "Trying XML User file : " << GetXMLUserFileName();

    if (xml_doc_user.ErrorID() < 1) {
      JSINFO << "Open XML User file : " << GetXMLUserFileName();
      xml_root_user =
          (tinyxml2::XMLElement *)xml_doc_user.FirstChildElement("jetscape");

      if (!xml_root_user) {
        JSWARN << "Not a valid JetScape XML User file!";
        exit(-1);
      }
    } else {
      JSWARN << "XML User file not found/not properly opened! Error code : "
             << xml_doc_user.ErrorID();
      exit(-1);
    }

    xml_user_file_open = true;
  }
}

//________________________________________________________________
void JetScapeXML::OpenXMLUserFile(string m_name) {
  SetXMLUserFileName(m_name);
  OpenXMLUserFile();
}

//________________________________________________________________
tinyxml2::XMLElement *
JetScapeXML::GetXMLElementMaster(std::initializer_list<const char *> &path) {

  JSWARN << "Deprecated function GetXMLElementMaster(). Call GetXMLElementMain() instead!";
  return GetXMLElementMain(path);
}

//________________________________________________________________
tinyxml2::XMLElement *
JetScapeXML::GetXMLElementMain(std::initializer_list<const char *> &path) {

  VERBOSE(2) << "Looking for element in Main file: " << path;

  OpenXMLMainFile();

  tinyxml2::XMLElement *currentElement = nullptr;
  for (auto &elementName : path) {
    if (!currentElement) {
      currentElement = xml_root_main->FirstChildElement(elementName);

      if (currentElement) {
        VERBOSE(3) << "Loaded " << elementName << " from xml_root_main";
      } else {
        VERBOSE(3) << elementName << " not found in xml_root_main";
      }
    } else {
      currentElement = currentElement->FirstChildElement(elementName);

      if (currentElement) {
        VERBOSE(3) << "Loaded " << elementName << " from "
                   << currentElement->Name();
      } else {
        VERBOSE(3) << elementName << " not found in " << elementName;
      }
    }
  }

  if (currentElement) {
    VERBOSE(2) << "Found element.";
  } else {
    VERBOSE(2) << "Did not find element.";
  }

  return currentElement;
}

//________________________________________________________________
tinyxml2::XMLElement *
JetScapeXML::GetXMLElementUser(std::initializer_list<const char *> &path) {

  VERBOSE(2) << "Looking for element in User file: " << path;

  OpenXMLUserFile();

  tinyxml2::XMLElement *currentElement = nullptr;
  for (auto &elementName : path) {
    if (!currentElement) {
      currentElement = xml_root_user->FirstChildElement(elementName);

      if (currentElement) {
        VERBOSE(3) << "Loaded " << elementName << " from xml_root_user";
      } else {
        VERBOSE(3) << elementName << " not found in xml_root_user";
      }
    } else {
      currentElement = currentElement->FirstChildElement(elementName);

      if (currentElement) {
        VERBOSE(3) << "Loaded " << elementName << " from "
                   << currentElement->Name();
      } else {
        VERBOSE(3) << elementName << " not found in " << elementName;
      }
    }
  }

  if (currentElement) {
    VERBOSE(2) << "Found element.";
  } else {
    VERBOSE(2) << "Did not find element.";
  }

  return currentElement;
}

//________________________________________________________________
tinyxml2::XMLElement *
JetScapeXML::GetElement(std::initializer_list<const char *> path,
                        bool isRequired /* = true */) {

  // Try to get value from User XML file
  tinyxml2::XMLElement *elementUser = GetXMLElementUser(path);
  if (elementUser) {
    return elementUser;
  }
  // Else, try to get value from Main XML file
  else {
    tinyxml2::XMLElement *elementMain = GetXMLElementMain(path);
    if (elementMain) {
      return elementMain;
    } else {
      if (isRequired) {
        JSWARN << "XML element " << path << " not found, but is required.";
        exit(-1);
      }
      return nullptr;
    }
  }
}

//________________________________________________________________
std::string
JetScapeXML::GetElementText(std::initializer_list<const char *> path,
                            bool isRequired /* = true */) {

  tinyxml2::XMLElement *element = GetElement(path, isRequired);

  if (element) {
    return element->GetText();
  } else {
    return "";
  }
}

//________________________________________________________________
int JetScapeXML::GetElementInt(std::initializer_list<const char *> path,
                               bool isRequired /* = true */) {

  tinyxml2::XMLElement *element = GetElement(path, isRequired);

  if (element) {
    int value = 0;
    element->QueryIntText(&value);
    return value;
  } else {
    return 0;
  }
}

//________________________________________________________________
double JetScapeXML::GetElementDouble(std::initializer_list<const char *> path,
                                     bool isRequired /* = true */) {

  tinyxml2::XMLElement *element = GetElement(path, isRequired);

  if (element) {
    double value = 0;
    element->QueryDoubleText(&value);
    return value;
  } else {
    return 0.;
  }
}

//________________________________________________________________
std::ostream &operator<<(std::ostream &os,
                         std::initializer_list<const char *> path) {

  int i = 0;
  int size = path.size();

  os << "\"";
  for (auto name : path) {
    os << name;
    if (i < size - 1) {
      os << ":";
    }
    i++;
  }
  os << "\"";

  return os;
}

} // end namespace Jetscape
