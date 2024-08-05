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

// Use CRTP for cloning of derived class in base class

#ifndef JETENERGYLOSSMODULE_H
#define JETENERGYLOSSMODULE_H

#include "JetEnergyLoss.h"

using std::abs;
using std::uniform_real_distribution;

namespace Jetscape {

template <typename Derived>
class JetEnergyLossModule : public JetEnergyLoss {
 public:
  using JetEnergyLoss::JetEnergyLoss;

  //! Deep copy constructor. The goal is to not have to initialize Modules
  //! multiple times
  virtual shared_ptr<JetEnergyLoss> Clone() const override {
    JSDEBUG << "Cloning task with id=" << GetId()
            << " and TaskNumber= " << GetMyTaskNumber();
    // DEBUG/TODO: KK: Joern's plan was to not have to call Init again, but I'm
    // not sure that can work/is desirable.
    auto ret = make_shared<Derived>(static_cast<const Derived &>(*this));
    // ret->Init();
    return ret;
  }

  //! Override deactivation
  void SetActive(bool m_active_exec) {
    throw std::runtime_error(
        "SetActive not supported for energy loss modules. "
        "Please remove the module from the manager.");
  };

 protected:
  /** Only one Eloss module at a time should be manipulating a parton
   * In the current setup, that's all but impossible to impose and relies
   * on cooperation between modules.
   * This is a crude way (relying on self-reporting) to check that this is
   * always the case.
   */
  bool TakeResponsibilityFor(Parton &p) {
    if (p.GetControlled()) {
      JSWARN << " Parton was controlled by " << p.GetController() << ". Now "
             << GetId() << " is trying to take responsibility as well.";
      throw std::runtime_error(
          "Two Eloss modules were fighting for one parton!");
    };
    // cout << " &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Was controlled by " <<
    // p.GetController() << endl;
    bool wascontrolled = p.SetController(GetId());
    // cout << " &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Now controlled by " <<
    // p.GetController() << endl;
    return wascontrolled;
  };
};

}  // end namespace Jetscape

#endif
