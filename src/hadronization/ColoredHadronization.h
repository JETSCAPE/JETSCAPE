#ifndef COLOREDHADRONIZATION_H
#define COLOREDHADRONIZATION_H

#include "HadronizationModule.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

class ColoredHadronization : public HadronizationModule<ColoredHadronization>
{  
 public:
  ColoredHadronization();
  virtual ~ColoredHadronization();
  
  void Init();
  void DoHadronization(vector<vector<shared_ptr<Parton>>>& shower, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);
    
private:
    double p_fake;
protected:
    static Pythia8::Pythia pythia;
};


#endif // COLOREDHADRONIZATION_H

