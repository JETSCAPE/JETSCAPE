// quick and dirty test class for Eloss modules ...

#ifndef ELOSSMODULESTEST_H
#define ELOSSMODULESTEST_H

#include "JetEnergyLossModule.h"

class Matter : public JetEnergyLossModule<Matter> //, std::enable_shared_from_this<Matter>
{  
 public:
  
  Matter();
  virtual ~Matter();

  void Init();
  void Exec();
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);
  
  //void SetQhat(double m_qhat) {qhat=m_qhat;}
  //double GetQhat() {return qhat;}
  
 private:

//double qhat;

};


class Martini : public JetEnergyLossModule<Martini> //, std::enable_shared_from_this<Martini>
{  
 public:
  
  Martini();
  virtual ~Martini();

  void Init();
  void Exec();
  virtual void WriteTask(weak_ptr<JetScapeWriter> w) {};
  
  //void SetQhat(double m_qhat) {qhat=m_qhat;}
  //double GetQhat() {return qhat;}
  
 private:

//double qhat;

};

#endif

