/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include <stdio.h>
#include <sys/stat.h>

#include <cstring>

#include "JetScapeLogger.h"
#include "FreestreamMilneWrapper.h"

using namespace std;

FreestreamMilneWrapper::FreestreamMilneWrapper() {
    preequilibrium_status_ = NOT_STARTED;
    SetId("Freestream-Milne");
}


FreestreamMilneWrapper::~FreestreamMilneWrapper() {
    if (preequilibrium_status_ != NOT_STARTED) delete fsmilne_ptr;
}


void FreestreamMilneWrapper::initialize_preequilibrium(PreEquilibriumParameterFile parameter_list) {
    INFO << "Initialize freestream-milne ...";
    VERBOSE(8);
    tinyxml2::XMLElement *para = GetPreequilibriumXML()->FirstChildElement("FreestreamMilne");
    if (!para) {
      WARN << " : freestream-milne not properly initialized in XML file ...";
      exit(-1);
    }
    string input_file = para->FirstChildElement("freestream_input_file")->GetText();//is this necessary? if we just force the user to have the 'freestream_input' file in the correct directory

    fsmilne_ptr = new FREESTREAMMILNE();
}


void FreestreamMilneWrapper::evolve_preequilibrium() {
    VERBOSE(8);
    INFO << "Initialize energy density profile in freestream-milne ...";
    // grab initial energy density from vector from initial state module
    std::vector<double> entropy_density = ini->entropy_density_distribution_; //note that this is the energy density when read by freestream-milne, not actually the entropy density!
    std::vector<float> entropy_density_float(entropy_density.begin(), entropy_density.end());
    fsmilne_ptr->initialize_from_vector(entropy_density_float);
    preequilibrium_status_ = INIT;
    if (preequilibrium_status_ == INIT) {
        INFO << "running freestream-milne ...";
        // evolve the medium via freestreaming
        fsmilne_ptr->run_freestream_milne();
        preequilibrium_status_ = DONE;
    }
    // now prepare to send the resulting hydro variables to the hydro module by coping hydro vectors to Preequilibrium base class members
    fsmilne_ptr->output_to_vectors(e_, P_, utau_, ux_, uy_, ueta_, pi00_, pi01_, pi02_, pi03_, pi11_, pi12_, pi13_, pi22_, pi23_, pi33_, bulk_Pi_);

}
