#ifndef MARTINIMUTEX_H
#define MARTINIMUTEX_H

#include <vector>
#include <string>
#include<iostream>
#include <memory>

#include "JetScapeTask.h"
#include "JetScapeModuleMutex.h"

using namespace Jetscape;
using std::shared_ptr;


class MartiniMutex : public JetScapeModuleMutex
{
  public:
    MartiniMutex();
    ~MartiniMutex();
    bool CheckMutex(vector<shared_ptr<JetScapeTask>> modules);  



};


#endif
