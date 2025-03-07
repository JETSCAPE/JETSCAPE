/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
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

#include <stdlib.h>

#include "JetScapeLogger.h"

using namespace std;

namespace Jetscape {

/**
 * @brief Singleton instance of the JetScapeXML class.
 */
JetScapeXML *JetScapeXML::m_pInstance = NULL;

/**
 * @brief Singleton instance accessor for JetScapeXML.
 *
 * This method returns the singleton instance of the JetScapeXML class.
 * If the instance does not exist, it creates a new one and logs the creation.
 *
 * @return JetScapeXML* Pointer to the singleton instance of JetScapeXML.
 */
JetScapeXML *JetScapeXML::Instance() {
  if (!m_pInstance) {
    JSINFO << "Created JetScapeXML Instance";
    m_pInstance = new JetScapeXML();
  }

  return m_pInstance;
}

/**
 * @brief Deprecated function to open the XML master file.
 *
 * This function is deprecated and should not be used. It logs a warning message
 * and calls the OpenXMLMainFile() function instead.
 *
 * @deprecated Use OpenXMLMainFile() instead.
 */
void JetScapeXML::OpenXMLMasterFile() {
  JSWARN << "Deprecated function OpenXMLMasterFile(): Call OpenXMLMainFile() "
            "instead!";
  OpenXMLMainFile();
}

/**
 * @brief Opens the main XML file for JetScape.
 *
 * This function attempts to open the main XML file specified by the
 * GetXMLMainFileName() function. If the file is successfully opened,
 * it sets the xml_root_main to the first "jetscape" element in the file.
 * If the main file cannot be opened, it attempts to open a deprecated
 * default master XML file located at "../config/jetscape_master.xml".
 * If neither file can be opened, the function logs a warning and exits
 * the program with an error code.
 *
 * @note The function sets the xml_main_file_open flag to true if the
 * main XML file is successfully opened.
 */
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
    } else {  // Check for an old default Master file
      auto errCode = xml_doc_main.ErrorID();  // Save the original error code
      SetXMLMainFileName(
          "../config/jetscape_master.xml");  // Try old default Master XML file
      xml_doc_main.LoadFile((char *)GetXMLMainFileName().c_str());
      VERBOSE(2) << "Looking for Old XML Master file : "
                 << GetXMLMainFileName();

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

/**
 * @brief Deprecated function to open the XML master file with a given name.
 *
 * This function is deprecated and should not be used. It logs a warning message
 * and calls the OpenXMLMainFile(string) function instead.
 *
 * @deprecated Use OpenXMLMainFile(string) instead.
 */
void JetScapeXML::OpenXMLMasterFile(string m_name) {
  JSWARN << "Deprecated function OpenXMLMasterFile(): Call OpenXMLMainFile() "
            "instead!";
  OpenXMLMainFile(m_name);
}

/**
 * @brief Opens the main XML file with the given name.
 *
 * This function sets the main XML file name using the provided string
 * and then opens the XML file.
 *
 * @param m_name The name of the XML file to be opened.
 */
void JetScapeXML::OpenXMLMainFile(string m_name) {
  SetXMLMainFileName(m_name);
  OpenXMLMainFile();
}

/**
 * @brief Opens the XML user file and loads its content.
 *
 * This function attempts to open and load the XML user file specified by
 * GetXMLUserFileName(). If the file is successfully opened and parsed, it
 * sets the xml_root_user to the first "jetscape" element in the file. If the
 * file cannot be opened or parsed, or if the root element is not found, the
 * function logs an appropriate warning message and terminates the program.
 *
 * The function ensures that the XML user file is only opened once by checking
 * the xml_user_file_open flag.
 *
 * @note This function will terminate the program if the XML user file cannot
 * be opened or parsed correctly, or if the root element is not found.
 */
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

/**
 * @brief Opens an XML user file.
 *
 * This function sets the XML user file name using the provided string
 * and then opens the XML user file.
 *
 * @param m_name The name of the XML user file to be opened.
 */
void JetScapeXML::OpenXMLUserFile(string m_name) {
  SetXMLUserFileName(m_name);
  OpenXMLUserFile();
}

/**
 * @brief Deprecated function to get the XML element from the specified path.
 *
 * This function is deprecated and should not be used. Use GetXMLElementMain()
 * instead.
 *
 * @deprecated Use GetXMLElementMain() instead.
 *
 * @param path An initializer list of const char* representing the path to the
 * XML element.
 * @return A pointer to the tinyxml2::XMLElement at the specified path.
 */
tinyxml2::XMLElement *JetScapeXML::GetXMLElementMaster(
    std::initializer_list<const char *> &path) {
  JSWARN << "Deprecated function GetXMLElementMaster(). Call "
            "GetXMLElementMain() instead!";
  return GetXMLElementMain(path);
}

/**
 * @brief Retrieves an XML element from the main XML file based on the provided
 * path.
 *
 * This function navigates through the XML structure starting from the root
 * element and follows the sequence of element names provided in the path to
 * locate the desired XML element.
 *
 * @param path A list of element names representing the path to the desired XML
 * element.
 * @return A pointer to the located tinyxml2::XMLElement, or nullptr if the
 * element is not found.
 */
tinyxml2::XMLElement *JetScapeXML::GetXMLElementMain(
    std::initializer_list<const char *> &path) {
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

/**
 * @brief Retrieves an XML element from the user file based on the provided
 * path.
 *
 * This function navigates through the XML structure using the given path and
 * returns the corresponding XML element if found. The path is specified as an
 * initializer list of element names.
 *
 * @param path An initializer list of const char* representing the path to the
 *             desired XML element.
 * @return A pointer to the tinyxml2::XMLElement if found, otherwise nullptr.
 */
tinyxml2::XMLElement *JetScapeXML::GetXMLElementUser(
    std::initializer_list<const char *> &path) {
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

/**
 * @brief Retrieves an XML element from the user or main XML file based on the
 * provided path.
 *
 * This function attempts to find an XML element specified by the given path. It
 * first searches in the user XML file, and if not found, it searches in the
 * main XML file. If the element is required and not found in either file, the
 * function logs a warning and terminates the program.
 *
 * @param path An initializer list of const char* representing the path to the
 * desired XML element.
 * @param isRequired A boolean flag indicating whether the XML element is
 * required. Defaults to true.
 * @return A pointer to the found tinyxml2::XMLElement, or nullptr if the
 * element is not found and isRequired is false.
 */
tinyxml2::XMLElement *JetScapeXML::GetElement(
    std::initializer_list<const char *> path, bool isRequired /* = true */) {
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

/**
 * @brief Retrieves the text content of an XML element specified by a path.
 *
 * This function navigates through the XML structure using the provided path
 * and returns the text content of the target element. If the element is not
 * found and isRequired is true, an exception may be thrown by GetElement.
 * If the element is not found and isRequired is false, an empty string is
 * returned.
 *
 * @param path An initializer list of const char* representing the path to the
 * target element.
 * @param isRequired A boolean flag indicating whether the element is required.
 * Defaults to true.
 * @return A std::string containing the text content of the target XML element,
 * or an empty string if the element is not found and isRequired is false.
 */
std::string JetScapeXML::GetElementText(
    std::initializer_list<const char *> path, bool isRequired /* = true */) {
  tinyxml2::XMLElement *element = GetElement(path, isRequired);

  if (element) {
    return element->GetText();
  } else {
    return "";
  }
}

/**
 * @brief Retrieves the int value of an XML element specified by the given
 * path.
 *
 * This function navigates through the XML structure using the provided path and
 * retrieves the int value of the target element.
 *
 * @param path An initializer list of const char pointers representing the path
 * to the target element.
 * @param isRequired A boolean flag indicating whether the element is required.
 * Defaults to true.
 * @return An int value representing the content of the target XML element, or
 * 0. if the element is not found.
 */
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

/**
 * @brief Retrieves the double value of an XML element specified by the given
 * path.
 *
 * This function navigates through the XML structure using the provided path and
 * retrieves the double value of the target element.
 *
 * @param path An initializer list of const char pointers representing the path
 * to the target element.
 * @param isRequired A boolean flag indicating whether the element is required.
 * Defaults to true.
 * @return A double value representing the content of the target XML element, or
 * 0. if the element is not found.
 */
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

/**
 * @brief Overloaded stream insertion operator for outputting an initializer
 * list of C-style strings. Prints the XML element path name.
 *
 * @param os The output stream to which the initializer list will be written.
 * @param path The initializer list of C-style strings to be output.
 * @return A reference to the output stream after the initializer list has been
 * written.
 */
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

}  // end namespace Jetscape
