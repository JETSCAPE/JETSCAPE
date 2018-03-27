#ifndef COLORLESSHADRONIZATION_H
#define COLORLESSHADRONIZATION_H

#include "HadronizationModule.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

class ColorlessHadronization : public HadronizationModule<ColorlessHadronization>
{  
 public:

  ColorlessHadronization();
  virtual ~ColorlessHadronization();
  
  void Init();
  void DoHadronization(vector<vector<shared_ptr<Parton>>>& shower, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);

 protected:
  static Pythia8::Pythia pythia;  
  
};


#endif // COLORLESSHADRONIZATION_H

