#ifndef ADSCFTMUTEX_H
#define ADSCFTMUTEX_H

#include <vector>
#include <string>
#include <iostream>
#include <memory>

#include "JetScapeTask.h"
#include "JetScapeModuleMutex.h"

using namespace Jetscape;
using std::shared_ptr;

class AdSCFTMutex : public JetScapeModuleMutex {
public:
  AdSCFTMutex();
  ~AdSCFTMutex();
  bool CheckMutex(vector<shared_ptr<JetScapeTask>> modules);
};

#endif
