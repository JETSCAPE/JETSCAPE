#ifndef CascadeTest_H
#define CascadeTest_H

#include "Afterburner.h"

using namespace Jetscape;

class CascadeTest : public Afterburner {
private:
  
public:

  CascadeTest() {SetId("CascadeTest");};
  virtual ~CascadeTest() {};

  void InitTask();

  virtual void CalculateTime();
  virtual void ExecTime();

  virtual void InitPerEvent();
  virtual void FinishPerEvent();

  virtual any GetHistory(); 

};

#endif // CascadeTest_H