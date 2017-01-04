// quick and dirty test class for Eloss modules ...

#ifndef ELOSSMODULESTEST_H
#define ELOSSMODULESTEST_H

#include "JetEnergyLoss.h"

class Matter : public JetEnergyLoss //, std::enable_shared_from_this<Matter>
{  
 public:
  
  Matter();
  virtual ~Matter();

  void Init();
  void Exec();

  //void SetQhat(double m_qhat) {qhat=m_qhat;}
  //double GetQhat() {return qhat;}
  
 private:

//double qhat;

};


class Martini : public JetEnergyLoss //, std::enable_shared_from_this<Martini>
{  
 public:
  
  Martini();
  virtual ~Martini();

  void Init();
  void Exec();
  
  //void SetQhat(double m_qhat) {qhat=m_qhat;}
  //double GetQhat() {return qhat;}
  
 private:

//double qhat;

};

#endif

