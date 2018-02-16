
#include <stdio.h>
#include <sys/stat.h>

#include <cstring>

#include "JetScapeLogger.h"
#include "freestream-milne_jetscape.h"

using namespace std;

FREESTREAM::FREESTREAM() {
  preequilibrium_status = NOT_START;
  SetId("Freestream-Milne");
}


FREESTREAM::~FREESTREAM() {
  if (preequilibrium_status != NOT_START) delete fsmilne_ptr;
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
  std::vector<double> entropy_density = ini->entropy_density_distribution_; //change this to the energy density!!!
  fsmilne_ptr->initialize_from_vector(entropy_density); //this needs to to be the energy density!!!
  preequilibrium_status = INITIALIZED;
  if (preequilibrium_status == INITIALIZED) {
    INFO << "running freestream-milne ...";
    //evolve the medium via freestreaming
    fsmilne_ptr->run_freestream_milne();
    preequilibrium_status = FINISHED;
  }

  //now send the resulting hydro variables to the hydro module

}

}
