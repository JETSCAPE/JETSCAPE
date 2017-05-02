// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: LongGang Pang (2017)
//                (UC Berkeley and LBNL)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include <boost/bind.hpp>
#include <boost/tokenizer.hpp>
#include "InitialCondition.h"

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
VarMap create_varmap(str projectile, str target,
                double cross_section, double grid_max, double grid_step)
{
  std::string cmd;
  cmd = projectile + " " + target
        + " -x " + std::to_string(cross_section)
        + " --grid-max " + std::to_string(grid_max)
        + " --grid-step " + std::to_string(grid_step);

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

void write_stream(std::ostream& os, int width,
    int num, double impact_param, const Event& event) {
  using std::fixed;
  using std::setprecision;
  using std::setw;
  using std::scientific;

  // Write a nicely-formatted line of event properties.
  os << setprecision(10)
     << setw(width)            << num
     << setw(15) << fixed      << impact_param
     << setw(5)                << event.npart()
     << setw(18) << scientific << event.multiplicity()
     << fixed;

  for (const auto& ecc : event.eccentricity())
    os << setw(14)             << ecc.second;

  for (const auto& phi_n : event.participant_plane())
    os << setw(14)             << phi_n.second;

  // Write the mass center (x, y)
  os << setw(14) << event.mass_center_index().first;
  os << setw(14) << event.mass_center_index().second;

  os << '\n';
}

} // end unnamedspace

// See header for explanation.
JetScapeInitial::~JetScapeInitial() = default;

// get one random collision in centrality for the given system
// stored_system = "auau200", "pbpb2760" or "pbpb5020"
// centrality = "0-5", "5-10", "30-40" or any range "a-b" for 0<=a<b<=100
JetScapeInitial::JetScapeInitial(str stored_system,
                    double centrality_low, double centrality_high,
                    double grid_max, double grid_step)
{
    VarMap var_map{};
    if (stored_system == "auau200") {
        str projectile = "Au";
        str target = "Au";
        double cross_section = 4.23;
        var_map = create_varmap(projectile, target, cross_section,
                grid_max, grid_step);
    } else if (stored_system == "pbpb2760") {
        str projectile = "Pb";
        str target = "Pb";
        double cross_section = 6.4;
        var_map = create_varmap(projectile, target, cross_section,
                grid_max, grid_step);
    } else if (stored_system == "pbpb5020") {
        str projectile = "Pb";
        str target = "Pb";
        double cross_section = 6.4;
        var_map = create_varmap(projectile, target, cross_section,
                grid_max, grid_step);
    }
    JetScapeCollision collision_(var_map);

    // smin, smax is stored as table
    collision_.sample_(100, 150);

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
}


JetScapeInitial::JetScapeInitial(str projectile, str target,
                double cross_section, double grid_max, double grid_step)
{
    VarMap var_map = create_varmap(projectile, target, cross_section,
                grid_max, grid_step);
    JetScapeCollision collision_(var_map);
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
}
