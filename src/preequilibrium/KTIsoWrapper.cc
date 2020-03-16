/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include <stdio.h>
#include <sys/stat.h>

#include <cstring>

#include "JetScapeLogger.h"
#include "KTIsoWrapper.h"

using namespace std;

KTIsoWrapper::KTIsoWrapper() {
    preequilibrium_status_ = NOT_STARTED;
    SetId("KTIso");
}


KTIsoWrapper::~KTIsoWrapper() {
    if (preequilibrium_status_ != NOT_STARTED) delete ktiso_ptr;
}


void KTIsoWrapper::InitializePreequilibrium(PreEquilibriumParameterFile parameter_list) {
    JSINFO << "Initialize KTIso ...";
    VERBOSE(8);
    tinyxml2::XMLElement *para = GetPreequilibriumXML()->FirstChildElement("FreestreamMilne");
    if (!para) {
      JSWARN << " : KTIso not properly initialized in XML file ...";
      exit(-1);
    }
    string input_file = para->FirstChildElement("freestream_input_file")->GetText();//is this necessary? if we just force the user to have the 'freestream_input' file in the correct directory

    ktiso_ptr = new ITA();
}


void KTIsoWrapper::EvolvePreequilibrium() {
    VERBOSE(8);
    JSINFO << "Initialize energy density profile in KTIso ...";
    // grab initial energy density from vector from initial state module
    std::vector<double> entropy_density = ini->GetEntropyDensityDistribution(); //note that this is the energy density when read by freestream-milne, not actually the entropy density!
    std::vector<float> entropy_density_float(entropy_density.begin(), entropy_density.end());
    ktiso_ptr->initialize_from_vector(entropy_density_float);
    preequilibrium_status_ = INIT;
    if (preequilibrium_status_ == INIT) {
        JSINFO << "running KTIso ...";
        // evolve the medium via freestreaming
        ktiso_ptr->run_ita();
        preequilibrium_status_ = DONE;
    }
    //set the value of end time of preequil dynamics for other modules to grab
    preequilibrium_tau_max_ = ktiso_ptr->tau_LandauMatch;
    JSINFO << "Setting end time of preequilibrium module to " << preequilibrium_tau_max_;
    // now prepare to send the resulting hydro variables to the hydro module by coping hydro vectors to Preequilibrium base class members
    ktiso_ptr->output_to_vectors(e_, P_, utau_, ux_, uy_, ueta_, pi00_, pi01_, pi02_, pi03_, pi11_, pi12_, pi13_, pi22_, pi23_, pi33_, bulk_Pi_);

}
