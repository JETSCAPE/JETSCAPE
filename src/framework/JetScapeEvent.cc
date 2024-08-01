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

#include "JetScapeEvent.h"

#include <iostream>

using namespace std;

namespace Jetscape {

JetScapeEvent::JetScapeEvent() {}

JetScapeEvent::JetScapeEvent(const JetScapeEvent &c) {
  partonCollection.clear();
  const vector<Parton> tmp = c.getPartonCollection();
  for (unsigned int ipart = 0; ipart < tmp.size(); ipart++) {
    partonCollection.push_back(c.getParton(ipart));
  }
}

JetScapeEvent::~JetScapeEvent() { partonCollection.clear(); }

const vector<Parton> &JetScapeEvent::getPartonCollection() const {
  return partonCollection;
}

const Parton &JetScapeEvent::getParton(int idx) const {
  return partonCollection.at(idx);
}

void JetScapeEvent::addParton(Parton &p) { partonCollection.push_back(p); }

void JetScapeEvent::addPartonShower(shared_ptr<PartonShower> ps) {
  for (unsigned int ipart = 0; ipart < ps->GetNumberOfPartons(); ipart++) {
    partonCollection.push_back(*(ps->GetPartonAt(ipart)));
  }
}

void JetScapeEvent::deleteParton(int idx) {
  partonCollection.erase(partonCollection.begin() +
                         idx);  // inefficient delete!!
}

}  // end namespace Jetscape
