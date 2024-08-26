#ifndef LBTMUTEX_H
#define LBTMUTEX_H

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "JetScapeModuleMutex.h"
#include "JetScapeTask.h"

using namespace Jetscape;
using std::shared_ptr;

class LBTMutex : public JetScapeModuleMutex {
 public:
  LBTMutex();
  ~LBTMutex();
  bool CheckMutex(vector<shared_ptr<JetScapeTask>> modules);
};

#endif
