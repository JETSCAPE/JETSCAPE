// quick and dirty test class for Eloss modules ...

#ifndef ELOSSMODULESTEST_H
#define ELOSSMODULESTEST_H

#include "JetEnergyLossModule.h"
#include "JetScapeWriter.h"
#include "JetScapeLogger.h"

#include <memory>

class Matter : public JetEnergyLossModule<Matter>  //, public std::enable_shared_from_this<Matter>
{  
 public:
  
  Matter();
  virtual ~Matter();

  void Init();
  void Exec();
  virtual void WriteTask(weak_ptr<JetScapeWriter> w) {DEBUG<<"TEST "<<GetId();};
  
  //void SetQhat(double m_qhat) {qhat=m_qhat;}
  //double GetQhat() {return qhat;}
  
 private:

//double qhat;

};


class Martini : public JetEnergyLossModule<Martini> // , public std::enable_shared_from_this<Martini>
{  
 public:
  
  Martini();
  virtual ~Martini();

  void Init();
  void Exec();
  void WriteTask(weak_ptr<JetScapeWriter> w) {};
  
  //void SetQhat(double m_qhat) {qhat=m_qhat;}
  //double GetQhat() {return qhat;}
  
 private:

//double qhat;

};

#endif

