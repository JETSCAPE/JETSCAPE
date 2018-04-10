/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * Modular, task-based framework
 * Intial Design: Joern Putschke, Kolja Kauder (Wayne State University)
 * For the full list of contributors see AUTHORS.

 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <boost/bind.hpp>
#include <boost/tokenizer.hpp>
#include "TrentoInitial.h"

namespace Jetscape {

// Helper functions for Collider ctor.
namespace {
/// @brief Tokenize a string.  The tokens will be separated by each non-quoted
///        space or equal character.  Empty tokens are removed.
///
/// @param input The string to tokenize.
///
/// @return Vector of tokens.
std::vector<std::string> tokenize(const std::string& input)
{
  typedef boost::escaped_list_separator<char> separator_type;
  separator_type separator("\\",    // The escape characters.
                           "= ",    // The separator characters.
                           "\"\'"); // The quote characters.

  // Tokenize the intput.
  boost::tokenizer<separator_type> tokens(input, separator);

  // Copy non-empty tokens from the tokenizer into the result.
  std::vector<std::string> result;
  copy_if(tokens.begin(), tokens.end(), std::back_inserter(result), 
          !boost::bind(&std::string::empty, _1));
  return result;
}

/// Auxiliary functions to create varmap from collision system
VarMap create_varmap(std::string projectile, std::string target,
                double cross_section, double grid_max, double grid_step,
                unsigned seed)
{
  VERBOSE(2) << "seed in create_varmap=" << seed;
  std::string cmd;
  cmd = projectile + " " + target
        + " -x " + std::to_string(cross_section)
        + " --grid-max " + std::to_string(grid_max)
        + " --grid-step " + std::to_string(grid_step)
        + " --random-seed " + std::to_string(seed);

  // Parse options with boost::program_options.
  // There are quite a few options, so let's separate them into logical groups.
  using OptDesc = po::options_description;

  using VecStr = std::vector<std::string>;
  OptDesc main_opts{};
  main_opts.add_options()
    ("projectile", po::value<VecStr>()->required()->
     notifier(  // use a lambda to verify there are exactly two projectiles
         [](const VecStr& projectiles) {
           if (projectiles.size() != 2)
            throw po::required_option{"projectile"};
           }),
     "projectile symbols")
    ("number-events", po::value<int>()->default_value(1),
     "number of events");

  // Make all main arguments positional.
  po::positional_options_description positional_opts{};
  positional_opts
    .add("projectile", 2)
    .add("number-events", 1);

  using VecPath = std::vector<fs::path>;
  OptDesc general_opts{"general options"};
  general_opts.add_options()
    ("help,h", "show this help message and exit")
    ("version", "print version information and exit")
    ("bibtex", "print bibtex entry and exit")
    // ("default-config", "print a config file with default settings and exit")
    ("config-file,c", po::value<VecPath>()->value_name("FILE"),
     "configuration file\n(can be passed multiple times)");

  OptDesc output_opts{"output options"};
  output_opts.add_options()
    ("quiet,q", po::bool_switch(),
     "do not print event properties to stdout")
    ("output,o", po::value<fs::path>()->value_name("PATH"),
     "HDF5 file or directory for text files");

  OptDesc phys_opts{"physical options"};
  phys_opts.add_options()
    ("reduced-thickness,p",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "reduced thickness parameter")
    ("fluctuation,k",
     po::value<double>()->value_name("FLOAT")->default_value(1., "1"),
     "gamma fluctuation shape parameter")
    ("nucleon-width,w",
     po::value<double>()->value_name("FLOAT")->default_value(.5, "0.5"),
     "Gaussian nucleon width [fm]")
    ("cross-section,x",
     po::value<double>()->value_name("FLOAT")->default_value(6.4, "6.4"),
     "inelastic nucleon-nucleon cross section sigma_NN [fm^2]")
    ("normalization,n",
     po::value<double>()->value_name("FLOAT")->default_value(1., "1"),
     "normalization factor")
    ("b-min",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "minimum impact parameter [fm]")
    ("b-max",
     po::value<double>()->value_name("FLOAT")->default_value(-1., "auto"),
     "maximum impact parameter [fm]")
    ("random-seed",
     po::value<int64_t>()->value_name("INT")->default_value(-1, "auto"),
     "random seed");

  OptDesc grid_opts{"grid options"};
  grid_opts.add_options()
    ("grid-max",
     po::value<double>()->value_name("FLOAT")->default_value(10., "10.0"),
     "xy max [fm]\n(grid extends from -max to +max)")
    ("grid-step",
     po::value<double>()->value_name("FLOAT")->default_value(0.2, "0.2"),
     "step size [fm]");

  // Make a meta-group containing all the option groups except the main
  // positional options (don't want the auto-generated usage info for those).
  OptDesc usage_opts{};
  usage_opts
    .add(general_opts)
    .add(output_opts)
    .add(phys_opts)
    .add(grid_opts);

  // Now a meta-group containing _all_ options.
  OptDesc all_opts{};
  all_opts
    .add(usage_opts)
    .add(main_opts);

  // Will be used several times.
  const std::string usage_str{
    "usage: trento [options] projectile projectile [number-events = 1]\n"};

  // Initialize a VarMap (boost::program_options::variables_map).
  // It will contain all configuration values.
  VarMap var_map{};

  // Parse command line options.
  //po::store(po::command_line_parser(argc, argv)
  //    .options(all_opts).positional(positional_opts).run(), var_map);

  po::store(po::command_line_parser(tokenize(cmd))
      .options(all_opts).positional(positional_opts).run(), var_map);

  // po::notify(var_map);
  return var_map;
}

//void write_stream(std::ostream& os, int width,
//    int num, double impact_param, const Event& event) {
//  using std::fixed;
//  using std::setprecision;
//  using std::setw;
//  using std::scientific;
//
//  // Write a nicely-formatted line of event properties.
//  os << setprecision(10)
//     << setw(width)            << num
//     << setw(15) << fixed      << impact_param
//     << setw(5)                << event.npart()
//     << setw(18) << scientific << event.multiplicity()
//     << fixed;
//
//  for (const auto& ecc : event.eccentricity())
//    os << setw(14)             << ecc.second;
//
//  for (const auto& phi_n : event.participant_plane())
//    os << setw(14)             << phi_n.second;
//
//  // Write the mass center (x, y)
//  os << setw(14) << event.mass_center_index().first;
//  os << setw(14) << event.mass_center_index().second;
//
//  os << '\n';
//}

} // end unnamedspace


// See header for explanation.
TrentoInitial::~TrentoInitial() = default;

// get one random collision in centrality for the given system
// stored_system = "AuAu200", "PbPb2760" or "PbPb5020"
// centrality = "0-5", "5-10", "30-40" or any range "a-b" for 0<=a<b<=100
void TrentoInitial::pre_defined(std::string stored_system,
                    double centrality_low, double centrality_high,
                    double grid_max, double grid_step, unsigned random_seed)
{
    VarMap var_map{};
    if (stored_system == "auau200") {
        std::string projectile = "Au";
        std::string target = "Au";
        double cross_section = 4.23;
        var_map = create_varmap(projectile, target, cross_section,
                grid_max, grid_step, random_seed);
    } else if (stored_system == "pbpb2760") {
        std::string projectile = "Pb";
        std::string target = "Pb";
        double cross_section = 6.4;
        var_map = create_varmap(projectile, target, cross_section,
                grid_max, grid_step, random_seed);
    } else if (stored_system == "pbpb5020") {
        std::string projectile = "Pb";
        std::string target = "Pb";
        double cross_section = 7.0;
        var_map = create_varmap(projectile, target, cross_section,
                grid_max, grid_step, random_seed);
    }

    TrentoCollision collision_(var_map);

    double smin, smax;
    std::tie(smin, smax) = get_entropy_range_(stored_system,
                                      centrality_low, centrality_high);
    collision_.sample_(smin, smax);
    info_ = collision_.info_;
    for (const auto& row : collision_.event_.reduced_thickness_grid()) {
      auto&& iter = row.begin();
      // Write all row elements except the last with a space delimiter afterwards.
      do {
        entropy_density_distribution_.push_back(*iter);
      } while (++iter != --row.end());
      // Write the last element and a linebreak.
      entropy_density_distribution_.push_back(*iter);
    }

    int nx = int(std::sqrt(entropy_density_distribution_.size()));
    double xmax = nx * grid_step / 2;
    SetRanges(xmax, xmax, 0.0);
    SetSteps(grid_step, grid_step, 0.0);
    compute_nbc();
}


void TrentoInitial::user_defined(std::string projectile, std::string target,
                double cross_section, double grid_max, double grid_step,
                unsigned random_seed)
{
    VarMap var_map = create_varmap(projectile, target, cross_section,
                grid_max, grid_step, random_seed);
    TrentoCollision collision_(var_map);
    double smin = 0; 
    double smax = std::numeric_limits<double>::max();
    collision_.sample_(smin, smax);
    info_ = collision_.info_;

    for (const auto& row : collision_.event_.reduced_thickness_grid()) {
      auto&& iter = row.begin();
      // Write all row elements except the last with a space delimiter afterwards.
      do {
        entropy_density_distribution_.push_back(*iter);
      } while (++iter != --row.end());
      // Write the last element and a linebreak.
      entropy_density_distribution_.push_back(*iter);
    }

    int nx = int(std::sqrt(entropy_density_distribution_.size()));
    VERBOSE(2) << "nx = " << nx;
    double xmax = nx * grid_step / 2;
    SetRanges(xmax, xmax, 0.0);
    SetSteps(grid_step, grid_step, 0.0);
    compute_nbc();
}


void TrentoInitial::InitTask() {
    INFO << " : Create initial condition ";
    trento_xml_ = xml_->FirstChildElement("Trento");
}


void TrentoInitial::Exec() {  
    VERBOSE(2) << " : Excute initial condition ";
    if (!trento_xml_) {
        WARN << " : Not a valid JetScape IS::Trento XML section in file!";
        exit(-1);
    } else {
        // trento_xml_->Attribute("A", "B") checks whether the attribute "A" has value "B"
        auto random_seed = (*get_mt19937_generator())();
        VERBOSE(2) << "Random seed used for TrentoInitial class" << random_seed;
        if ( trento_xml_->Attribute("use_module", "pre_defined") ) {
            auto predef = trento_xml_->FirstChildElement("pre_defined");
            std::string collision_system(predef->Attribute("collision_system"));
            VERBOSE(2) << "collision_system=" << collision_system;
            double centrality_min = std::atof(predef->Attribute("centrality_min"));
            double centrality_max = std::atof(predef->Attribute("centrality_max"));
            pre_defined(collision_system, centrality_min, centrality_max,
                    get_x_max(), get_x_step(), random_seed);
        } else if (trento_xml_->Attribute("use_module", "user_defined") ) {
            auto usrdef = trento_xml_->FirstChildElement("user_defined");
            std::string projectile(usrdef->Attribute("projectile"));
            std::string target(usrdef->Attribute("target"));
            // center of mass collision energy per pair of nucleon
            double sqrts_NN = std::atof(usrdef->Attribute("sqrts"));
            double cross_section = std::atof(usrdef->Attribute("cross_section"));
            user_defined(projectile, target, cross_section,
                         get_x_max(), get_x_step(), random_seed);
        }
    }
}

void TrentoInitial::Clear() {
    VERBOSE(2) << " : Finish creating initial condition ";
    entropy_density_distribution_.clear();
    num_of_binary_collisions_.clear();
}


TrentoInitial::TrentoInitial() : InitialState() {
    SetId("Trento");
}

/** Notice that this function assumes the total number of charged
 * particles is propotional to initial total entropy, this is true
 * for constant eta/s. this function returns the entropy range for
 * one given centrality class from reading a stored collision
 * system*/
std::tuple<double, double> TrentoInitial::get_entropy_range_(std::string collision_system,
        double centrality_low, double centrality_high) {
    std::stringstream centrality_class_path;
    centrality_class_path << "data_table/trento_" << collision_system
                          << "_cent.csv";

    std::ifstream fin(centrality_class_path.str());
    if (!fin.is_open()) {
        throw runtime_error("open "
                + centrality_class_path.str() + " failed");
    }
    char buf[256];
    fin.getline(buf, 256);

    double centrality_bound[101];
    double entropy_bound[101];

    for (int i = 0; i < 101; i++) {
        fin >> centrality_bound[i] >> entropy_bound[i];
    }

    fin.close();

    auto interp1d = [&](double cent) {
        if (cent < 0.0 || cent > 100.0) {
            throw RangeFailure(std::string("centrality ") + std::to_string(cent)
                    + " is not in the range [0, 100]");
        } else {
            double x0, x1, y0, y1;
            int j = int(cent);
            if (j > cent) {
                x0 = centrality_bound[j-1];
                y0 = entropy_bound[j-1];
                x1 = centrality_bound[j];
                y1 = entropy_bound[j];
            } else {
                x0 = centrality_bound[j];
                y0 = entropy_bound[j];
                x1 = centrality_bound[j+1];
                y1 = entropy_bound[j+1];
            }
            double r  = (cent - x0)/(x1 - x0);
            return (1.0 - r) * y0 + r * y1;
        }
    };

    return std::make_pair(interp1d(centrality_high),
                          interp1d(centrality_low));
}

void TrentoInitial::compute_nbc() {
    for ( const auto si : entropy_density_distribution_ ) {
        auto si_squre = si * si;
        // this works for IP-Glasma like initial condition
        num_of_binary_collisions_.push_back(si_squre);
    }
}



} // end namespace Jetscape
