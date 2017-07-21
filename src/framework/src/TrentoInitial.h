// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: LongGang Pang (2017)
//                (UC Berkeley and LBNL)
// Use part code from trento (because of the private property)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef TRENTO_INITIAL_H
#define TRENTO_INITIAL_H

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
#include "InitialState.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"

using OptDesc = po::options_description;

using namespace trento;

namespace Jetscape {

////////////////////////// Trento Initial Condition Wrapper //////////////////////
class TrentoInitial : public InitialState {
  public:
    // Initialize from XML configuration
    TrentoInitial();

    // get one random collision in centrality range 0-100%
    void user_defined(std::string projectile, std::string target,
                    double cross_section, double grid_max,
                    double grid_step);

    // get one random collision in centrality for the given system
    // stored_system = "auau200", "pbpb2760" or "pbpb5020"
    // centrality_range = [centrality_min, centrality_max]
    void pre_defined(std::string stored_system,
                    double centrality_min, double centrality_max,
                    double grid_max, double grid_step);

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

    std::tuple<double, double> get_entropy_range_(std::string collision_system,
        double centrality_low, double centrality_high);

    // compute number of binary collisions
    void compute_nbc();

    tinyxml2::XMLElement * trento_xml_;

    /// The output instance.
    // Output output_;
};

} // end namespace Jetscape

#endif
