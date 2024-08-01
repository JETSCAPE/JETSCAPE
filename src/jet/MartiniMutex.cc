#include "MartiniMutex.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "JetScapeTask.h"

using namespace std;
using std::shared_ptr;
using namespace Jetscape;

MartiniMutex::MartiniMutex() {}

MartiniMutex::~MartiniMutex() {}

bool MartiniMutex::CheckMutex(vector<shared_ptr<JetScapeTask>> modules) {
  bool isLbt = false;
  bool isAdscft = false;

  for (auto module : modules) {
    string name = module->GetId();
    if (!name.compare("LBT"))
      isLbt = true;
    if (!name.compare("AdSCFT"))
      isAdscft = true;
  }

  if (isLbt || isAdscft)
    return false;
  return true;
}
