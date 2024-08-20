#ifndef JETSCAPEMODULEMUTEX_H
#define JETSCAPEMODULEMUTEX_H

#include <memory>
#include <vector>

#include "JetScapeTask.h"

using namespace std;
using std::shared_ptr;

namespace Jetscape {

class JetScapeModuleMutex {
 public:
  virtual bool CheckMutex(vector<shared_ptr<JetScapeTask>> modules) = 0;
};

}  // end namespace Jetscape

#endif
