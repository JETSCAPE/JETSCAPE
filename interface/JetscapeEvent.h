/*******************************************************************************

This is the base data structure class of Jetscape, which will get passed to all
modules to read and modify as needed. This base class should evolve as new needs
arsise. 

*******************************************************************************/




#ifndef JetscapeEvent_h
#define JetscapeEvent_h

/*  Could use HepMC as well if we decide to include it  */
// #include "HepMC/GenEvent.h"
// #include "HepMC/GenParticle.h"
#include "JetClass.h"
#include "PartonClass.h"
#include "HydroClass.h"
#include "ConfigClass.h"

class JetScapeEvent {
 public :
  std::vector<Jet> jetCollection; //defined in JetClass.h
  std::vector<Parton> partonCollection; //defined in JetClass.h
  Hydro hydroCollection; //defined in HydroClass.h
  Config parameterSet;
}
#endif 
