/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/amajumder/JETSCAPE-COMP/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#ifndef JETSCAPEEVENT_H
#define JETSCAPEEVENT_H

#include "JetClass.h"
#include "PartonShower.h"

namespace Jetscape {

class JetScapeEvent
{
  
 public:
  
  JetScapeEvent();
  JetScapeEvent(const JetScapeEvent &c); //copy constructor
  ~JetScapeEvent();

  const Parton& getParton(int idx) const;
  const vector<Parton>& getPartonCollection() const;
  void addParton(Parton &p);
  void addPartonShower(shared_ptr<PartonShower> ps);
  void deleteParton(int idx);

 private:

  vector<Parton> partonCollection;
  
};

} // end namespace Jetscape

#endif
