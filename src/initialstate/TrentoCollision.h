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
// Use part code from trento (because of the private property)

#ifndef TRENTOCOLLISION_H
#define TRENTOCOLLISION_H

#include <tuple>
#include <memory>
#include <iostream>
#include <iomanip>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>

#include "fwd_decl.h"
#include "event.h"
#include "nucleon.h"
#include "output.h"

// CMake sets this definition.
// Fall back to a sane default.
#ifndef TRENTO_VERSION_STRING
#define TRENTO_VERSION_STRING "dev"
#endif

#include "JetScapeModuleBase.h"
#include "tinyxml2.h"

using str = std::string;

using OptDesc = po::options_description;

using namespace trento;

namespace Jetscape {

typedef struct {
    double impact_parameter;
    double num_participant;
    double total_entropy;
    std::map<int, double> ecc; // order, eccentricity
    std::map<int, double> psi; // order, participant_plane
    double xmid, ymid;
} EventInfo;

////////////////////////// Trento Initial Condition Wrapper //////////////////////
class TrentoCollision{
  public:
    //// get one random collision in centrality range 0-100%
    //TrentoCollision(str projectile, str target,
    //                double cross_section, double grid_max,
    //                double grid_step);

    //// get one random collision in centrality for the given system
    //// stored_system = "auau200", "pbpb2760" or "pbpb5020"
    //// centrality = "0-5", "5-10", "30-40" or any range "a-b"
    //TrentoCollision(str stored_system, str centrality,
    //                double grid_max, double grid_step);

    explicit TrentoCollision(const VarMap & var_map);

    /// Declare a destructor to properly handle the std::unique_ptr<Nucleus>
    ~TrentoCollision();

    // sample jet production position from Ta * Tb
    // where Ta * Tb is the distribution of num_of_binary collisions
    // return (x, y) tuple
    std::tuple<double, double> jet_production_position_();

    // return (bimp, Event) tuple whose entropy in (smin, smax)
    void sample_(double smin, double smax);

    /// The event instance.
    Event event_;

    EventInfo info_;

  private:
    // return var_map which holds options to run trento
    VarMap trento_options_(str cmd);

    // 2D entropy density distribution from sqrt(Ta * Tb)
    // std::vector<float> entropy_density_;

    /// Sample a min-bias impact parameter within the set range.
    double sample_impact_param();

    /// Pair of nucleus projectiles.
    std::unique_ptr<Nucleus> nucleusA_, nucleusB_;

    /// The nucleon profile instance.
    NucleonProfile nucleon_profile_;

    /// Number of events to run.
    const int nevents_;

    /// Minimum and maximum impact parameter.
    const double bmin_, bmax_;

    /// Parameterizes the degree of asymmetry between the two projectiles.  Used
    /// to apportion the total impact parameter to each projectile so that the
    /// resulting overlap is approximately centered.  Given the two nuclear radii
    /// rA and rB (see Nucleus::radius()), the parameter is
    ///
    ///   asymmetry = rA / (rA + rB)
    ///
    /// and the impact parameter b is apportioned
    ///
    ///   offset_A = asymmetry * b
    ///   offset_B = (asymmetry-1) * b
    ///
    /// So for example, two identical projectiles would have asymmetry = 0.5, and
    /// each would be offset by half the impact parameter.  A proton-lead system
    /// would have asymmetry = 1, meaning the lead nucleus would be offset by the
    /// entire impact parameter and the proton would not be offset at all.
    const double asymmetry_;


    /// The output instance.
    Output output_;
};

} // end namespace Jetscape

#endif // TRENTOCOLLISION_H
