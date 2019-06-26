// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// TRENTO3D: Three-dimensional extension of TRENTO by Weiyao Ke
// MIT License

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <climits>

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
	"TRENTO:\n"
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
    "}\n\n";

  std::cout <<
	"TRENTO3D:\n"
    "@article{Ke:2016jrd,\n"
    "      author         = \"Ke, Weiyao and Moreland, J. Scott and Bernhard, \n"
	"                         Jonah E. and Bass, Steffen A.\",\n"
    "      title          = \"{Constraints on rapidity-dependent initial conditions\n"
    "                        from charged particle pseudorapidity densities and\n"
    "                        two-particle correlations}\",\n"
    "      journal        = \"Phys. Rev.\",\n"
    "      volume         = \"C96\",\n"
    "      year           = \"2017\",\n"
    "      number         = \"4\",\n"
    "      pages          = \"044912\",\n"
    "      doi            = \"10.1103/PhysRevC.96.044912\",\n"
    "      eprint         = \"1610.08490\",\n"
    "      archivePrefix  = \"arXiv\",\n"
    "      primaryClass   = \"nucl-th\",\n"
    "      SLACcitation   = \"%%CITATION = ARXIV:1610.08490;%%\"\n";
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
     "HDF5 file or directory for text files")
    ("no-header", po::bool_switch(),
     "do not write headers to text files");

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
    ("nucleon-min-dist,d",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "minimum nucleon-nucleon distance [fm]")
    ("mean-coeff,m",
     po::value<double>()->value_name("FLOAT")->default_value(1., "1."),
     "rapidity mean coefficient")
    ("std-coeff,s",
     po::value<double>()->value_name("FLOAT")->default_value(3., "3."),
     "rapidity std coefficient")
    ("skew-coeff,t",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0."),
     "rapidity skew coefficient")
	 ("skew-type,r",
      po::value<int>()->value_name("INT")->default_value(1, "1"),
      "rapidity skew type: 1: relative, 2: absolute, other: no skew")
    ("jacobian,j",
     po::value<double>()->value_name("FLOAT")->default_value(0.8, "0.8"),
     "<pt>/<mt> used in Jacobian")
    ("normalization,n",
     po::value<double>()->value_name("FLOAT")->default_value(1., "1"),
     "normalization factor");

  OptDesc coll_opts{"collision options"};
  coll_opts.add_options()
    ("beam-energy,e",
     po::value<double>()->value_name("FLOAT")->default_value(2760, "2760"),
     "collision beam energy sqrt(s) [GeV], initializes cross section")
    ("cross-section,x",
     po::value<double>()->value_name("FLOAT")->default_value(-1, "off"),
     "manual inelastic nucleon-nucleon cross section sigma_NN [fm^2]")
    ("b-min",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "minimum impact parameter [fm]")
    ("b-max",
     po::value<double>()->value_name("FLOAT")->default_value(-1., "auto"),
     "maximum impact parameter [fm]")
	("npart-min",
     po::value<int>()->value_name("INT")->default_value(0, "0"),
     "minimum Npart cut")
    ("npart-max",
     po::value<int>()->value_name("INT")->default_value(
     std::numeric_limits<int>::max(), "INT_MAX"), "maximum Npart cut")
    ("s-min",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "minimum entropy cut")
    ("s-max",
     po::value<double>()->value_name("FLOAT")->default_value(
     std::numeric_limits<double>::max(), "DOUBLE_MAX"), "maxmimum entropy cut")
    ("random-seed",
     po::value<int64_t>()->value_name("INT")->default_value(-1, "auto"),
     "random seed")
    ("ncoll,b", po::bool_switch(),
     "calculate # of binary collision and binary collision density");

  OptDesc grid_opts{"grid options"};
  grid_opts.add_options()
    ("xy-max",
     po::value<double>()->value_name("FLOAT")->default_value(10., "10.0"),
     "xy max [fm]\n(transverse grid from -max to +max)")
    ("xy-step",
     po::value<double>()->value_name("FLOAT")->default_value(0.2, "0.2"),
     "transverse step size [fm]")
    ("eta-max",
     po::value<double>()->value_name("FLOAT")->default_value(0.0, "0.0"),
     "pseudorapidity max \n(eta grid from -max to +max)")
    ("eta-step",
     po::value<double>()->value_name("FLOAT")->default_value(0.5, "0.5"),
     "pseudorapidity step size");

  // Make a meta-group containing all the option groups except the main
  // positional options (don't want the auto-generated usage info for those).
  OptDesc usage_opts{};
  usage_opts
    .add(general_opts)
    .add(output_opts)
    .add(phys_opts)
    .add(coll_opts)
    .add(grid_opts);

  // Now a meta-group containing _all_ options.
  OptDesc all_opts{};
  all_opts
    .add(usage_opts)
    .add(main_opts);

  // Will be used several times.
  const std::string usage_str{
    "usage: trento [options] projectile projectile [number-events = 1]\n"};
  const std::string usage_str3d{
    "To operate in 3D mode, make sure --eta-max is nonzero.\n"};

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
        << usage_str << usage_str3d
        << "\n"
           "projectile = { p | d | Cu | Cu2 | Xe | Au | Au2 | Pb | U | U2 | U3 }\n"
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
  }
  catch (const po::required_option&) {
    // Handle this exception separately from others.
    // This occurs e.g. when the program is excuted with no arguments.
    std::cerr << usage_str << usage_str3d
			  << "run 'trento --help' for more information\n";
    return 1;
  }
  catch (const std::exception& e) {
    // For all other exceptions just output the error message.
    std::cerr << e.what() << '\n';
    return 1;
  }

  return 0;
}
