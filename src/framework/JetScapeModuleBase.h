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

#ifndef JETSCAPEMODULEBASE_H
#define JETSCAPEMODULEBASE_H

#include <map>
#include <memory>
#include <random>
#include <string>

#include "JetScapeTask.h"
#include "JetScapeXML.h"
#include "sigslot.h"

namespace Jetscape {

class JetScapeWriter;

class JetScapeModuleBase
    : public JetScapeTask,
      public sigslot::has_slots<sigslot::multi_threaded_local> {
 public:
  /** Default constructor to create a JetScapeModuleBase. It sets the XML file
   * name to a default string value.
   */
  JetScapeModuleBase();

  /** This is a constructor to create a JetScapeModuleBase. It sets the XML file
   * name to "m_name" to be used to read input parameters.
   */
  JetScapeModuleBase(string m_name);

  /** This is a destructor for the JetScapeModuleBase.
   */
  virtual ~JetScapeModuleBase();

  // virtual shared_ptr<JetScapeModuleBase> Clone() const {return nullptr;}
  /** A virtual function for a default initialization of a JetScapeModuleBase.
   * It also checks whether a XML file is loaded or not.
   */
  virtual void Init();

  /** A virtual function to define a default Exec() function for a
   * JetScapeModuleBase. It can be overridden by different modules/tasks.
   */
  virtual void Exec(){};

  /**  A virtual function to define a default Clear() function for a
   * JetScapeModuleBase. It can be overridden by different modules/tasks.
   */
  virtual void Clear(){};

  /** This function sets the name of the XML file to be used to store output
   * information for the modules/tasks of a JetScapeTask.
   */
  void SetXMLMainFileName(string m_name) { xml_main_file_name = m_name; }

  /** This function returns the XML file name. This file contains the output
   * data for the modules/tasks of a JetScapeTask.
   */
  string GetXMLMainFileName() { return xml_main_file_name; }

  /** This function sets the name of the XML file to be used to store output
   * information for the modules/tasks of a JetScapeTask.
   */
  void SetXMLUserFileName(string m_name) { xml_user_file_name = m_name; }

  /** This function returns the XML file name. This file contains the output
   * data for the modules/tasks of a JetScapeTask.
   */
  string GetXMLUserFileName() { return xml_user_file_name; }

  /** This function returns the current event number.
   */
  static int GetCurrentEvent() { return current_event; }

  /** This function increases the current event number by one.
   */
  static void IncrementCurrentEvent() { current_event++; }

  /** This function returns a random number based on Mersenne-Twister algorithm.
   */
  shared_ptr<std::mt19937> GetMt19937Generator();

  /** Helper functions for XML parsing, wrapping functionality in JetScapeXML:
   */
  tinyxml2::XMLElement *GetXMLElement(std::initializer_list<const char *> path,
                                      bool isRequired = true) {
    return JetScapeXML::Instance()->GetElement(path, isRequired);
  }
  std::string GetXMLElementText(std::initializer_list<const char *> path,
                                bool isRequired = true) {
    return JetScapeXML::Instance()->GetElementText(path, isRequired);
  }
  int GetXMLElementInt(std::initializer_list<const char *> path,
                       bool isRequired = true) {
    return JetScapeXML::Instance()->GetElementInt(path, isRequired);
  }
  double GetXMLElementDouble(std::initializer_list<const char *> path,
                             bool isRequired = true) {
    return JetScapeXML::Instance()->GetElementDouble(path, isRequired);
  }

 private:
  std::string xml_main_file_name;
  std::string xml_user_file_name;
  static int current_event;
  shared_ptr<std::mt19937> mt19937_generator_;
};

/**
 * @class JetScapeModuleComponentFactory
 * @brief Factory for modules in the Jetscape framework.
 *
 * This class implements a static map (i.e. shared between all instances of
 * JetScapeModuleBase) consisting of a std::string of module names (i.e. name in
 * XML config) and a function that creates an instance of the module. This will
 * allow us to automatically add new modules to Jetscape without modifying the
 * framework classes.
 *
 * Based on: https://stackoverflow.com/a/582456 and
 * https://github.com/alisw/AliPhysics/blob/master/PWG/EMCAL/EMCALtasks/AliEmcalCorrectionComponent.h
 */

/// Template function for creating a new module. Used to register the module.
template <typename T>
shared_ptr<JetScapeModuleBase> createT() {
  return std::make_shared<T>();
}

// Factory to create and keep track of new modules
class JetScapeModuleFactory {
 public:
  virtual ~JetScapeModuleFactory() {}

  typedef std::map<std::string, shared_ptr<JetScapeModuleBase> (*)()> map_type;

  /// Creates an instance of an object based on the name if the name is
  /// registered in the map.
  static shared_ptr<JetScapeModuleBase> createInstance(std::string const &s) {
    map_type::iterator it = getMap()->find(s);
    if (it == getMap()->end()) {
      return 0;
    }
    return it->second();
  }

 protected:
  /// Creates and access the module map
  static map_type *getMap() {
    // We never delete the map (until program termination) because we cannot
    // guarantee correct destruction order
    if (!moduleMap) {
      moduleMap = new map_type;
    }
    return moduleMap;
  }

 private:
  /// Contains the map to all of the modules
  static map_type *moduleMap;
};

/**
 * @class RegisterJetScapeModule
 * @brief Registers Jetscape modules in the factory map
 */
template <typename T>
class RegisterJetScapeModule : public JetScapeModuleFactory {
 public:
  /// Registers the name of the module to map to a function that can create the
  /// module
  RegisterJetScapeModule(std::string const &s) {
    getMap()->insert(std::make_pair(s, &createT<T>));
  }
};

}  // end namespace Jetscape

#endif
