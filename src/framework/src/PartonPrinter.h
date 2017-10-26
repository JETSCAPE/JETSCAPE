#ifndef PARTONPRINTER_H
#define PARTONPRINTER_H

#include "JetScapeModuleBase.h"
#include "PartonShower.h"
#include <vector>
#include <string>


namespace Jetscape {

class PartonPrinter : public JetScapeModuleBase
{

public:

PartonPrinter();
virtual ~PartonPrinter();

virtual void Init();
virtual void Exec() final;
virtual void Clear();

void GetFinalPartons(shared_ptr<PartonShower> pShower, vector<shared_ptr<Parton>>& fPartons);

void GetFinalPartons2(shared_ptr<PartonShower> pShower);

void GetPartonsAtTime(shared_ptr<PartonShower> pShower, vector<shared_ptr<Parton>>& fPartons, double time);

void PrintFinalPartons(vector<shared_ptr<Parton>>& fPartons) {fPartons = pFinals;};

private:

vector<shared_ptr<Parton>> pFinals;

};

} // end namespace Jetscape

#endif
