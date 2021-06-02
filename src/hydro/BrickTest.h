#ifndef BrickTest_H
#define BrickTest_H

#include "Brick.h"

using namespace Jetscape;

class BrickTest : public Brick {
private:
  
public:

  BrickTest() : Brick() {SetId("BrickTest");};
  virtual ~BrickTest() {};

  virtual void CalculateTime();
  virtual void ExecTime();

  virtual void InitPerEvent();
  virtual void FinishPerEvent();

  virtual any GetHistory(); 

};

#endif // BrickTest_H