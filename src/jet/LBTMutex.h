#ifndef LBTMUTEX_H
#define LBTMUTEX_H

#include <vector>
#include <string>
#include<iostream>
#include <memory>

#include "JetScapeTask.h"
#include "JetScapeModuleMutex.h"

using namespace Jetscape;
using std::shared_ptr;


class LBTMutex : public JetScapeModuleMutex
{
  public:
    LBTMutex();
    ~LBTMutex();
    bool CheckMutex(vector<shared_ptr<JetScapeTask>> modules);



};


#endif
