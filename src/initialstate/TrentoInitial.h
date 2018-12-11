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

#ifndef TRENTOINITIAL_H
#define TRENTOINITIAL_H

#include <tuple>
#include <memory>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>

#include "fwd_decl.h"
#include "JetScapeModuleBase.h"
#include "tinyxml2.h"
#include "collider.h"
#include "InitialState.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"

using OptDesc = po::options_description;
using VarMap = po::variables_map;
using namespace trento;

namespace Jetscape {

typedef struct {
    double impact_parameter;
    double num_participant;
    double num_binary_collisions;
    double total_entropy;
    std::map<int, double> ecc; // order, eccentricity
    std::map<int, double> psi; // order, participant_plane
    double xmid, ymid;
} EventInfo;

/**The output data format (from http://qcd.phy.duke.edu/trento/usage.html#output-options):
 * The grid will always be a square N × N array, with N = ceil(2*max/step).
 * So e.g. the default settings (max = 10 fm, step = 0.2 fm) imply a 100 × 100 grid.
 * The ceiling function ensures that the number of steps is always rounded up,
 * so e.g. given max = 10 fm and step 0.3 fm, the grid will be 67 × 67.
 * In this case, the actual grid max will be marginally increased (max = nsteps*step/2).
**/

////////////////////////// Trento Initial Condition Wrapper //////////////////////
class TrentoInitial : public InitialState {
  public:
    // Initialize from XML configuration
    TrentoInitial();
    ~TrentoInitial();

    //void Init();
    void Exec();
    void Clear();
    void InitTask();

    struct RangeFailure : public std::runtime_error {
        using std::runtime_error::runtime_error;
    };
	EventInfo info_;

  private:

    tinyxml2::XMLElement * trento_xml_;
	std::shared_ptr<trento::Collider> TrentoGen_;
	std::pair<double, double> GenCenTab(std::string proj, std::string targ, VarMap var_map, int cL, int cH);
    /// The output instance.
    // Output output_;
};

} // end namespace Jetscape

#endif // TRENTOINITIAL_H
