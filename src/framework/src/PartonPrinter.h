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

vector<shared_ptr<Parton>> GetFinalPartons(shared_ptr<PartonShower> pShower);

vector<shared_ptr<Parton>> GetPartonsAtTime(shared_ptr<PartonShower> pShower, double time);

void PrintFinalPartons(shared_ptr<PartonShower> pShower);

//shared_ptr<PartonShower> pShower;
//vector<shared_ptr<Parton>> partons;


//static shared_ptr<PartonShower> getShower(){return pShower;}
//static void setShower(shared_ptr<PartonShower> shower){pShower = shower;}

//static void setPartons(vector<shared_ptr<Parton>> vPartons){partons = vPartons;}
//static vector<shared_ptr<Parton>> getPartons(){return  partons;}

//shared_ptr<PartonShower> pShower;

vector<shared_ptr<Parton>> partons;


};

} // end namespace Jetscape

#endif
