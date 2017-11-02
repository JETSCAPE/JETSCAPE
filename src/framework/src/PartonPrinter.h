#ifndef PARTONPRINTER_H
#define PARTONPRINTER_H

#include "JetScapeModuleBase.h"
#include "PartonShower.h"
#include <vector>
#include <string>
#include<fstream>


namespace Jetscape {

class PartonPrinter : public JetScapeModuleBase
{

public:

PartonPrinter();
virtual ~PartonPrinter();

virtual void Init();
virtual void Exec() final;
virtual void Clear();
    std::ofstream dist_output; ///< the output stream where events are saved to file

void GetFinalPartons(shared_ptr<PartonShower> pShower, vector<shared_ptr<Parton>>& fPartons);

void GetFinalPartons2(shared_ptr<PartonShower> pShower);

void GetPartonsAtTime(shared_ptr<PartonShower> pShower, vector<shared_ptr<Parton>>& fPartons, double time);

void PrintFinalPartons(vector<shared_ptr<Parton>>& fPartons) {fPartons = pFinals;};

private:

vector<shared_ptr<Parton>> pFinals;

};

} // end namespace Jetscape

#endif
