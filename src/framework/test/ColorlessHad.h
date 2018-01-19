#ifndef COLORLESSHAD_H
#define COLORLESSHAD_H

#include "HadronizationModule.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

class ColorlessHad : public HadronizationModule<ColorlessHad>
{  
 public:

  ColorlessHad();
  virtual ~ColorlessHad();
  
  void Init();
  void DoHadronization(vector<vector<shared_ptr<Parton>>>& shower, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);

};


#endif



















