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

#ifndef JETSCAPEEVENT_H
#define JETSCAPEEVENT_H

#include "JetClass.h"
#include "PartonShower.h"

namespace Jetscape {

/**
 * @class JetScapeEvent 
 * @brief Class for events in JETSCAPE.
 * 
 * The event class contains a collection of partons and parton showers.
 */
class JetScapeEvent {
 public:
  /**
   * @brief Constructor for JetScapeEvent.
   */
  JetScapeEvent();

  /**
  * @brief Copy constructor for JetScapeEvent.
  * 
  * @param c The JetScapeEvent object to copy.
  * 
  * Copies the parton collection from another JetScapeEvent instance.
  */
  JetScapeEvent(const JetScapeEvent &c);  // copy constructor
  
  /**
  * @brief Destructor for JetScapeEvent.
  */
  ~JetScapeEvent();

  /**
  * @brief Get a parton from the collection.
  * 
  * @param idx The index of the parton to retrieve.
  * @return The parton at the specified index.
  */
  const Parton &getParton(int idx) const;
  
  /**
  * @brief Get the parton collection.
  *
  * @return The parton collection.
  */
  const vector<Parton> &getPartonCollection() const;

  /**
  * @brief Add a parton to the collection.
  * 
  * @param p The parton to add.
  */
  void addParton(Parton &p);

  /**
  * @brief Add a parton shower to the collection.
  * 
  * @param ps The parton shower to add.
  */
  void addPartonShower(shared_ptr<PartonShower> ps);

  /**
  * @brief Delete a parton from the collection.
  * 
  * @param idx The index of the parton to delete.
  * 
  * @todo Innefficient delete. TODO
  */
  void deleteParton(int idx);

 private:
  vector<Parton> partonCollection;
};

}  // end namespace Jetscape

#endif
