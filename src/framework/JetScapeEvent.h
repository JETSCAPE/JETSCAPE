// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef JETSCAPEEVENT_H
#define JETSCAPEEVENT_H

#include "JetClass.hpp"
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
