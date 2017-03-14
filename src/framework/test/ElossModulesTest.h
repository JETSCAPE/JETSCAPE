// quick and dirty test class for Eloss modules ...

#ifndef ELOSSMODULESTEST_H
#define ELOSSMODULESTEST_H

#include "JetEnergyLossModule.h"

class Matter : public JetEnergyLossModule<Matter> //, public std::enable_shared_from_this<Matter>
{  
 public:
  
  Matter();
  virtual ~Matter();

  void Init();
  void Exec();
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);
  
 private:

};


class Martini : public JetEnergyLossModule<Martini> //, public std::enable_shared_from_this<Martini>
{  
 public:
  
  Martini();
  virtual ~Martini();

  void Init();
  void Exec();
  virtual void WriteTask(weak_ptr<JetScapeWriter> w) {};
  
 private:

};

#endif

