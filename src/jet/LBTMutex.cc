#include <vector>
#include <string>
#include <iostream>
#include <memory>

#include "JetScapeTask.h"
#include "LBTMutex.h"

using namespace std;
using std::shared_ptr;
using namespace Jetscape;

LBTMutex::LBTMutex() {}

LBTMutex::~LBTMutex() {}

bool LBTMutex::CheckMutex(vector<shared_ptr<JetScapeTask>> modules) {
  bool isMartini = false;
  bool isAdscft = false;

  for (auto module : modules) {
    string name = module->GetId();
    if (!name.compare("Martini"))
      isMartini = true;
    if (!name.compare("AdSCFT"))
      isAdscft = true;
  }

  if (isMartini || isAdscft)
    return false;
  return true;
}
