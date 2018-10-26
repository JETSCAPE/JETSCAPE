#include <vector>
#include <string>
#include<iostream>
#include <memory>

#include "JetScapeTask.h"
#include "AdSCFTMutex.h"

using namespace std;
using std::shared_ptr;
using namespace Jetscape;

AdSCFTMutex::AdSCFTMutex()
{
}

AdSCFTMutex::~AdSCFTMutex()
{
}

bool AdSCFTMutex::CheckMutex(vector<shared_ptr<JetScapeTask>> modules)
{
  cout << "##### inside check mutex method of martini" << endl;
  bool isLbt = false;
  bool isMartini = false;

  for(auto module : modules)
  {
    string name = module->GetId();
    if(!name.compare("LBT")) isLbt = true;
    if(!name.compare("Martini")) isMartini = true;
  }

  if(isLbt || isMartini) return false;
  return true;
}
