
#ifndef TEST_FREESTREAM_JETSCAPE_H_
#define TEST_FREESTREAM_JETSCAPE_H_

#include "PreequilibriumDynamics.h"
#include "FreestreamMilne.cpp" //this file is in 3rdparty/freestream-milne/src/wrapper - how do we include it besides a hard link?

using namespace Jetscape;

class FREESTREAM: public PreequilibriumDynamics {
 private:
    int mode; //!< records running mode
    FREESTREAMMILNE *fsmilne_ptr;

 public:
    FREESTREAM();
    ~FREESTREAM();

    void initialize_preequilibrium(PreEquilibriumParameterFile parameter_list);

    void evolve_preequilibrium();
};

//JUST FOR TESTING! REMOVE THIS BLOCK ONCE CMAKE LINKING FILES WORKS !
/*
FREESTREAM::FREESTREAM() {
  preequilibrium_status = NOT_STARTED;
  SetId("Freestream-Milne");
}
FREESTREAM::~FREESTREAM() {
  if (preequilibrium_status != NOT_STARTED) delete fsmilne_ptr;
}

void FREESTREAM::initialize_preequilibrium(Parameter parameter_list) {
  INFO << "Initialize freestream-milne ...";
  VERBOSE(8);
  tinyxml2::XMLElement *para = GetPreequilibriumXML()->FirstChildElement("Preequilibrium");
  if (!para) {
    WARN << " : freestream-milne not properly initialized in XML file ...";
    exit(-1);
  }
  string input_file = para->FirstChildElement("freestream_milne_input_file")->GetText(); //is this necessary? if we just force the user to have the 'freestream_input' file in the correct directory

  fsmilne_ptr = new FREESTREAMMILNE();

}

void FREESTREAM::evolve_preequilibrium() {
  VERBOSE(8);
  INFO << "Initialize energy density profile in freestream-milne ...";
  //grab initial energy density from vector from initial state module
  std::vector<double> entropy_density = ini->entropy_density_distribution_; //note that this is the energy density when read by freestream-milne
  std::vector<float> entropy_density_float(entropy_density.begin(), entropy_density.end()); //converting to a float for now
  fsmilne_ptr->initialize_from_vector(entropy_density_float);
  preequilibrium_status = INIT;
  if (preequilibrium_status == INIT) {
    INFO << "running freestream-milne ...";
    //evolve the medium via freestreaming
    fsmilne_ptr->run_freestream_milne();
    preequilibrium_status = DONE;
  }
  //now send the resulting hydro variables to the hydro module
}
*/
//JUST FOR TESTING! REMOVE THIS BLOCK ONCE CMAKE LINKING FILES WORKS !

#endif  // TEST_FREESTREAM_JETSCAPE_H_
