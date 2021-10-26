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

#ifndef HADRONIZATIONMODULE_H
#define HADRONIZATIONMODULE_H

#include "Hadronization.h"

using std::uniform_real_distribution;

namespace Jetscape {

template <typename Derived> class HadronizationModule : public Hadronization {

public:
  using Hadronization::Hadronization;

  virtual shared_ptr<Hadronization> Clone() const override {
    auto ret = make_shared<Derived>(static_cast<const Derived &>(*this));
    return ret;
  }
};

} // end namespace Jetscape

#endif
