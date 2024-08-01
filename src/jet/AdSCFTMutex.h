#ifndef ADSCFTMUTEX_H
#define ADSCFTMUTEX_H

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "JetScapeModuleMutex.h"
#include "JetScapeTask.h"

using namespace Jetscape;
using std::shared_ptr;

class AdSCFTMutex : public JetScapeModuleMutex {
 public:
  AdSCFTMutex();
  ~AdSCFTMutex();
  bool CheckMutex(vector<shared_ptr<JetScapeTask>> modules);
};

#endif
