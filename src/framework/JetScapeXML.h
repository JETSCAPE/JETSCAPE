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

#ifndef JETSCAPEXML_H
#define JETSCAPEXML_H

#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <string>

#include "tinyxml2.h"

using std::runtime_error;
using std::string;

namespace Jetscape {

/**
 * @class JetScapeXML
 * @brief JetScape XML init reader class (meant as singleton)
 *
 * This class contains the machinery to load two XML configuration files: a Main
 * file, and a User file.
 *
 */
class JetScapeXML {
 public:
  static JetScapeXML *Instance();

  /**
   * @note Master file: These functions are deprecated. Users should use
   * the Main functions instead. These functions have been updated to use
   * the 'main' instead of 'master' variables
   */

  /**
   * @brief Get the root element of the master XML file.
   * @return The root element of the master XML file.
   * @deprecated Use GetXMLRootMain() instead.
   */
  tinyxml2::XMLElement *GetXMLRootMaster() { return xml_root_main; }

  /**
   * @brief Get the XML document of the master XML file.
   * @return The XML document of the master XML file.
   * @deprecated Use GetXMLDocumentMain() instead.
   */
  tinyxml2::XMLDocument &GetXMLDocumentMaster() { return xml_doc_main; }

  /**
   * @brief Get an XML element from the master XML file.
   * @param path The path to the XML element.
   * @return The XML element at the specified path.
   * @deprecated Use GetXMLElementMain() instead.
   */
  tinyxml2::XMLElement *GetXMLElementMaster(
      std::initializer_list<const char *> &path);

  /**
   * @brief Set the name of the master XML file.
   * @param m_name The name of the master XML file.
   * @deprecated Use SetXMLMainFileName() instead.
   */
  void SetXMLMasterFileName(string m_name) { xml_main_file_name = m_name; }

  /**
   * @brief Get the name of the master XML file.
   * @return The name of the master XML file.
   * @deprecated Use GetXMLMainFileName() instead.
   */
  std::string GetXMLMasterFileName() { return xml_main_file_name; }

  /**
   * @brief Check if the master XML file is open.
   * @return True if the master XML file is open, false otherwise.
   * @deprecated Use IsMainFileOpen() instead.
   */
  bool IsMasterFileOpen() { return xml_main_file_open; }

  /**
   * @brief Open the master XML file.
   * @deprecated Use OpenXMLMainFile() instead.
   */
  void OpenXMLMasterFile();

  /**
   * @brief Open the master XML file with the specified name.
   * @param m_name The name of the master XML file.
   * @deprecated Use OpenXMLMainFile() instead.
   */
  void OpenXMLMasterFile(string m_name);

  // Main file:

  /**
   * @brief Get the root element of the main XML file.
   * @return The root element of the main XML file.
   */
  tinyxml2::XMLElement *GetXMLRootMain() const { return xml_root_main; }

  /**
   * @brief Get the XML document of the main XML file.
   * @return The XML document of the main XML file.
   */
  tinyxml2::XMLDocument &GetXMLDocumentMain() { return xml_doc_main; }

  /**
   * @brief Get an XML element from the main XML file.
   * @param path The path to the XML element.
   * @return The XML element at the specified path.
   */
  tinyxml2::XMLElement *GetXMLElementMain(
      std::initializer_list<const char *> &path);

  /**
   * @brief Set the name of the main XML file.
   * @param m_name The name of the main XML file.
   */
  void SetXMLMainFileName(string m_name) { xml_main_file_name = m_name; }

  /**
   * @brief Get the name of the main XML file.
   * @return The name of the main XML file.
   */
  std::string GetXMLMainFileName() const { return xml_main_file_name; }

  /**
   * @brief Check if the main XML file is open.
   * @return True if the main XML file is open, false otherwise.
   */
  bool IsMainFileOpen() const { return xml_main_file_open; }

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
  void OpenXMLMainFile();

  /**
   * @brief Opens the main XML file with the given name.
   *
   * This function sets the main XML file name using the provided string
   * and then opens the XML file.
   *
   * @param m_name The name of the XML file to be opened.
   */
  void OpenXMLMainFile(string m_name);

  // User file

  /**
   * @brief Get the root element of the user XML file.
   * @return The root element of the user XML file.
   */
  tinyxml2::XMLElement *GetXMLRootUser() { return xml_root_user; }

  /**
   * @brief Get the XML document of the user XML file.
   * @return The XML document of the user XML file.
   */
  tinyxml2::XMLDocument &GetXMLDocumentUser() { return xml_doc_user; }

  /**
   * @brief Get an XML element from the user XML file.
   * @param path The path to the XML element.
   * @return The XML element at the specified path.
   */
  tinyxml2::XMLElement *GetXMLElementUser(
      std::initializer_list<const char *> &path);

  /**
   * @brief Set the name of the user XML file.
   * @param m_name The name of the user XML file.
   */
  void SetXMLUserFileName(string m_name) { xml_user_file_name = m_name; }

  /**
   * @brief Get the name of the user XML file.
   * @return The name of the user XML file.
   */
  std::string GetXMLUserFileName() { return xml_user_file_name; }

  /**
   * @brief Check if the user XML file is open.
   * @return True if the user XML file is open, false otherwise.
   */
  bool IsUserFileOpen() { return xml_user_file_open; }

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
  void OpenXMLUserFile();

  /**
   * @brief Opens an XML user file.
   *
   * This function sets the XML user file name using the provided string
   * and then opens the XML user file.
   *
   * @param m_name The name of the XML user file to be opened.
   */
  void OpenXMLUserFile(string m_name);

  /**
   * @note Helper functions for XML parsing
   * Look first in user XML file for a parameter, and if not found
   * look in the main XML file.
   */

  /**
   * @brief Retrieves an XML element from the user or main XML file based on the
   * provided path.
   *
   * This function attempts to find an XML element specified by the given path.
   * It first searches in the user XML file, and if not found, it searches in
   * the main XML file. If the element is required and not found in either file,
   * the function logs a warning and terminates the program.
   *
   * @param path An initializer list of const char* representing the path to the
   * desired XML element.
   * @param isRequired A boolean flag indicating whether the XML element is
   * required. Defaults to true.
   * @return A pointer to the found tinyxml2::XMLElement, or nullptr if the
   * element is not found and isRequired is false.
   */
  tinyxml2::XMLElement *GetElement(std::initializer_list<const char *> path,
                                   bool isRequired = true);

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
   * @param isRequired A boolean flag indicating whether the element is
   * required. Defaults to true.
   * @return A std::string containing the text content of the target XML
   * element, or an empty string if the element is not found and isRequired is
   * false.
   */
  std::string GetElementText(std::initializer_list<const char *> path,
                             bool isRequired = true);

  /**
   * @brief Retrieves the int value of an XML element specified by the given
   * path.
   *
   * This function navigates through the XML structure using the provided path
   * and retrieves the int value of the target element.
   *
   * @param path An initializer list of const char pointers representing the
   * path to the target element.
   * @param isRequired A boolean flag indicating whether the element is
   * required. Defaults to true.
   * @return An int value representing the content of the target XML element, or
   * 0. if the element is not found.
   */
  int GetElementInt(std::initializer_list<const char *> path,
                    bool isRequired = true);

  /**
   * @brief Retrieves the double value of an XML element specified by the given
   * path.
   *
   * This function navigates through the XML structure using the provided path
   * and retrieves the double value of the target element.
   *
   * @param path An initializer list of const char pointers representing the
   * path to the target element.
   * @param isRequired A boolean flag indicating whether the element is
   * required. Defaults to true.
   * @return A double value representing the content of the target XML element,
   * or 0. if the element is not found.
   */
  double GetElementDouble(std::initializer_list<const char *> path,
                          bool isRequired = true);

 private:
  /**
   * @brief Default constructor for JetScapeXML.
   *
   * This constructor initializes the JetScapeXML class with default values.
   */
  JetScapeXML() {
    xml_main_file_name = "";
    xml_main_file_open = false;
    xml_user_file_name = "";
    xml_user_file_open = false;
  };

  /**
   * @brief Copy constructor for JetScapeXML.
   *
   * This constructor creates a new instance of JetScapeXML by copying an
   * existing instance.
   *
   * @param other The JetScapeXML instance to copy from.
   */
  JetScapeXML(JetScapeXML const &){};

  /**
   * @brief Singleton instance of the JetScapeXML class.
   */
  static JetScapeXML *m_pInstance;

  // Main file

  /**
   * @note use unique pointer here instead of raw pointer
   * (check with tinyxml interface)
   */
  tinyxml2::XMLElement *xml_root_main;
  tinyxml2::XMLDocument xml_doc_main;

  std::string xml_main_file_name;
  bool xml_main_file_open;

  // User file

  /**
   * @note use unique pointer here instead of raw pointer
   * (check with tinyxml interface)
   */
  tinyxml2::XMLElement *xml_root_user;
  tinyxml2::XMLDocument xml_doc_user;

  std::string xml_user_file_name;
  bool xml_user_file_open;
};

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
                         std::initializer_list<const char *> path);

}  // end namespace Jetscape

#endif
