// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>

#include "collider.h"
#include "fwd_decl.h"

// CMake sets this definition.
// Fall back to a sane default.
#ifndef TRENTO_VERSION_STRING
#define TRENTO_VERSION_STRING "dev"
#endif

namespace trento {

namespace {

void print_version() {
  std::cout << "trento " << TRENTO_VERSION_STRING << '\n';
}

void print_bibtex() {
  std::cout <<
    "@article{Moreland:2014oya,\n"
    "      author         = \"Moreland, J. Scott and Bernhard, Jonah E. and Bass,\n"
    "                        Steffen A.\",\n"
    "      title          = \"{Alternative ansatz to wounded nucleon and binary\n"
    "                        collision scaling in high-energy nuclear collisions}\",\n"
    "      journal        = \"Phys.Rev.\",\n"
    "      number         = \"1\",\n"
    "      volume         = \"C92\",\n"
    "      pages          = \"011901\",\n"
    "      doi            = \"10.1103/PhysRevC.92.011901\",\n"
    "      year           = \"2015\",\n"
    "      eprint         = \"1412.4708\",\n"
    "      archivePrefix  = \"arXiv\",\n"
    "      primaryClass   = \"nucl-th\",\n"
    "      SLACcitation   = \"%%CITATION = ARXIV:1412.4708;%%\",\n"
    "}\n";
}

// TODO
// void print_default_config() {
//   std::cout << "to do\n";
// }

}  // unnamed namespace

}  // namespace trento

int main(int argc, char* argv[]) {
  using namespace trento;

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

  try {
    // Initialize a VarMap (boost::program_options::variables_map).
    // It will contain all configuration values.
    VarMap var_map{};

    // Parse command line options.
    po::store(po::command_line_parser(argc, argv)
        .options(all_opts).positional(positional_opts).run(), var_map);

    // Handle options that imply immediate exit.
    // Must do this _before_ po::notify() since that can throw exceptions.
    if (var_map.count("help")) {
      std::cout
        << usage_str
        << "\n"
           "projectile = { p | d | Cu | Cu2 | Au | Au2 | Pb | U | U2 | U3 }\n"
        << usage_opts
        << "\n"
           "see the online documentation for complete usage information\n";
      return 0;
    }
    if (var_map.count("version")) {
      print_version();
      return 0;
    }
    if (var_map.count("bibtex")) {
      print_bibtex();
      return 0;
    }
    // if (var_map.count("default-config")) {
    //   print_default_config();
    //   return 0;
    // }

    // Merge any config files.
    if (var_map.count("config-file")) {
      // Everything except general_opts.
      OptDesc config_file_opts{};
      config_file_opts
        .add(main_opts)
        .add(output_opts)
        .add(phys_opts)
        .add(grid_opts);

      for (const auto& path : var_map["config-file"].as<VecPath>()) {
        if (!fs::exists(path)) {
          throw po::error{
            "configuration file '" + path.string() + "' not found"};
        }
        fs::ifstream ifs{path};
        po::store(po::parse_config_file(ifs, config_file_opts), var_map);
      }
    }

    // Save all the final values into var_map.
    // Exceptions may occur here.
    po::notify(var_map);

    // Go!
    Collider collider{var_map};
    collider.run_events();
  } catch (const po::required_option&) {
    // Handle this exception separately from others.
    // This occurs e.g. when the program is excuted with no arguments.
    std::cerr << usage_str << "run 'trento --help' for more information\n";
    return 1;
  } catch (const std::exception& e) {
    // For all other exceptions just output the error message.
    std::cerr << e.what() << '\n';
    return 1;
  }

  return 0;
}
