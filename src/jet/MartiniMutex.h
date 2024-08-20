#ifndef MARTINIMUTEX_H
#define MARTINIMUTEX_H

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "JetScapeModuleMutex.h"
#include "JetScapeTask.h"

using namespace Jetscape;
using std::shared_ptr;

class MartiniMutex : public JetScapeModuleMutex {
 public:
  MartiniMutex();
  ~MartiniMutex();
  bool CheckMutex(vector<shared_ptr<JetScapeTask>> modules);
};

#endif
