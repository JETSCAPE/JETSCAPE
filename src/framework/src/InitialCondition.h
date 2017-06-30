// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: LongGang Pang (2017)
//                (UC Berkeley and LBNL)
// Use part code from trento (because of the private property)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef INITIAL_CONDITION_H
#define INITIAL_CONDITION_H

#include <tuple>
#include <memory>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>

#include "fwd_decl.h"
#include "JetScapeModuleBase.h"
#include "tinyxml2.h"
#include "Collision.h"

using str = std::string;

using OptDesc = po::options_description;

using namespace trento;

namespace Jetscape {

////////////////////////// Trento Initial Condition Wrapper //////////////////////
class JetScapeInitial{
  public:
    // get one random collision in centrality range 0-100%
    JetScapeInitial(str projectile, str target,
                    double cross_section, double grid_max,
                    double grid_step);

    // get one random collision in centrality for the given system
    // stored_system = "auau200", "pbpb2760" or "pbpb5020"
    // centrality = "0-5", "5-10", "30-40" or any range "a-b"
    JetScapeInitial(str stored_system, double centrality_low, double centrality_high,
                    double grid_max, double grid_step);

    ~JetScapeInitial();

    // sample jet production position from Ta * Tb
    // where Ta * Tb is the distribution of num_of_binary collisions
    // return (x, y) tuple
    std::tuple<double, double> jet_production_position_();

    std::vector<double> entropy_density_distribution_;

    EventInfo info_;

    struct RangeFailure : public std::runtime_error {
        using std::runtime_error::runtime_error;
    };

  private:

    std::tuple<double, double> get_entropy_range_(str collision_system,
        double centrality_low, double centrality_high);

    /// The output instance.
    // Output output_;
};

} // end namespace Jetscape

#endif
