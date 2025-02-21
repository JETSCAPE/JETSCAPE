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

#ifndef HADRONIZATIONMODULE_H
#define HADRONIZATIONMODULE_H

#include "Hadronization.h"

namespace Jetscape {

/**
 * @class HadronizationModule
 * @brief Template class for hadronization modules in JETSCAPE.
 *
 * This class serves as a template for hadronization modules, allowing derived
 * classes to implement their specific hadronization models while supporting
 * dynamic cloning.
 *
 * @tparam Derived The specific hadronization model class that inherits from this template.
 */
template <typename Derived>
class HadronizationModule : public Hadronization {
 public:
  /// Inherit the constructor from the base class.
  using Hadronization::Hadronization;

  /**
   * @brief Creates a clone of the current hadronization module.
   *
   * This function dynamically allocates a new instance of the derived class
   * using the copy constructor and returns a shared pointer to it.
   *
   * @return A shared pointer to the cloned hadronization module.
   */
  virtual shared_ptr<Hadronization> Clone() const override {
    auto ret = make_shared<Derived>(static_cast<const Derived &>(*this));
    return ret;
  }
};

}  // end namespace Jetscape

#endif
