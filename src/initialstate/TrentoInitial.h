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
#include "TrentoCollision.h"
#include "InitialState.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"

using OptDesc = po::options_description;

using namespace trento;

namespace Jetscape {

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

    // get one random collision in centrality range 0-100%
    void UserDefined(std::string projectile, std::string target,
                    double cross_section, double normalisation, double grid_max,
                    double grid_step, unsigned random_seed);

    // get one random collision in centrality for the given system
    // stored_system = "auau200", "pbpb2760" or "pbpb5020"
    // centrality_range = [centrality_min, centrality_max]
    void PreDefined(std::string stored_system, double normalisation,
                    double centrality_min, double centrality_max,
                    double grid_max, double grid_step, unsigned random_seed);

    ~TrentoInitial();

    //void Init();
    void Exec();
    void Clear();
    void InitTask();

    EventInfo info_;

    struct RangeFailure : public std::runtime_error {
        using std::runtime_error::runtime_error;
    };

  private:

    std::tuple<double, double> GetEntropyRange(std::string collision_system,
        double centrality_low, double centrality_high);

    // compute number of binary collisions
    void ComputeNbc();

    tinyxml2::XMLElement * trento_xml_;

    /// The output instance.
    // Output output_;
};

} // end namespace Jetscape

#endif // TRENTOINITIAL_H
