#ifndef Pythia8_h
#define Pythia8_h

#include <"BaseHardScattering.">
#include <Pythia8/Pythia.h>
#include <Pythia8Plugins/HepMC2.h>


class Pythia8 : BaseHardScattering {

  public :
    Pythia8(Parameter parameter_list);
    ~Pythia8();
    
    void generatePartons(HepMC::GenEvent & event) override;

}

#endif
