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

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <boost/bind.hpp>
#include <boost/tokenizer.hpp>
#include <algorithm>
#include <functional>
#include <string>
#include <sstream>
#include "JetScapeLogger.h"

#include "TrentoInitial.h"

namespace Jetscape {

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
}// end namespace

// See header for explanation.
TrentoInitial::TrentoInitial() : InitialState() {
    SetId("Trento");
}

TrentoInitial::~TrentoInitial() = default;

void TrentoInitial::InitTask() {
    JSINFO << " Initialzie TRENTo initial condition ";
    trento_xml_ = xml_->FirstChildElement("Trento");


    if (!trento_xml_) {
        JSWARN << " : Not a valid JetScape IS::Trento XML section in file!";
        exit(-1);
    }
	// TRENTO OPTION DESK
	using namespace trento;
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

	// NOW LETS FILL IN THE OPTION DESK
	auto phy_opts = trento_xml_->FirstChildElement("PhysicsInputs");
	auto cut_opts = trento_xml_->FirstChildElement("CutInputs");
	auto trans_opts = trento_xml_->FirstChildElement("TransInputs");
	auto longi_opts = trento_xml_->FirstChildElement("LongiInputs");

	double xymax = GetXMax(), dxy = GetXStep();
	double etamax = GetZMax(), deta = GetZStep();
	
	auto random_seed = (*GetMt19937Generator())();
	//TEMPORARY FOR TESTING
	//auto random_seed = 1;
	//TEMPORARY
	JSINFO << "Random seed used for Trento " << random_seed;

	std::string proj(phy_opts->Attribute("projectile"));
	std::string targ(phy_opts->Attribute("target"));
	double sqrts = std::atof(phy_opts->Attribute("sqrts"));
	double cross_section = std::atof(phy_opts->Attribute("cross-section"));
	double normalization = std::atof(phy_opts->Attribute("normalization"));
	
	int cen_low = std::atoi(cut_opts->Attribute("centrality-low"));
	int cen_high = std::atoi(cut_opts->Attribute("centrality-high"));

	double p = std::atof(trans_opts->Attribute("reduced-thickness"));
	double k = std::atof(trans_opts->Attribute("fluctuation"));
	double w = std::atof(trans_opts->Attribute("nucleon-width"));
	double d = std::atof(trans_opts->Attribute("nucleon-min-dist"));


	double mean = std::atof(longi_opts->Attribute("mean-coeff"));
	double var = std::atof(longi_opts->Attribute("std-coeff"));
	double skew = std::atof(longi_opts->Attribute("skew-coeff"));
	int skew_type = std::atof(longi_opts->Attribute("skew-type"));
	double J = std::atof(longi_opts->Attribute("jacobian"));
        
        std::ostringstream stream1;
        stream1 << std::fixed  
		<< " --random-seed " << random_seed
		<< " --cross-section " << cross_section
		<< " --beam-energy " << sqrts
		<< " --reduced-thickness " << p
		<< " --fluctuation " << k
		<< " --nucleon-width " << w
		<< " --nucleon-min-dist " << d
		<< " --mean-coeff " << mean
		<< " --std-coeff " << var 
		<< " --skew-coeff " << skew
		<< " --skew-type " << skew_type
		<< " --jacobian " << J
		<< " --quiet ";
	std::string options1 = stream1.str();
        std::ostringstream stream2;
	stream2 << std::fixed	
                << " --normalization " << normalization
		<< " --ncoll " // calcualte # of binary collision
		<< std::setprecision(20) << " --xy-max " << xymax
		<< std::setprecision(20) << " --xy-step " << dxy
		<< std::setprecision(20) << " --eta-max " << etamax
		<< std::setprecision(20) << " --eta-step " << deta;
        std::string options2 = stream2.str();
	// Handle centrality table, not normzlized, default grid, 2D (fast) !!!
	std::string cmd_basic = proj+" "+targ+" 10000 "+options1;
	VarMap var_map_basic{}; 
  	po::store(po::command_line_parser(tokenize(cmd_basic))
      .options(all_opts).positional(positional_opts).run(), var_map_basic);
   
        std::string options_cut = "";    
        if (cen_low == 0 && cen_high == 100) {
             JSINFO << "TRENTo Minimum Biased Mode Generates 0-100(%) of nuclear inelastic cross-section";
        }
        else {
      	    auto Ecut = GenCenTab(proj, targ, var_map_basic, cen_low, cen_high);
	    double Ehigh = Ecut.first*normalization; // rescale the cut
	    double Elow = Ecut.second*normalization; // rescale the cut
        
	    JSINFO << "The total energy density cut for centrality = [" << cen_low << ", "
	   	 << cen_high << "] (%) is:";
	    JSINFO << Elow << "<dE/deta(eta=0)<" << Ehigh;
	    options_cut = 
		  " --s-max " + std::to_string(Ehigh)
		+ " --s-min " + std::to_string(Elow);
	    // Set trento configuration
	}
	std::string cmd = proj+" "+targ+" 1 "+options1+options2+options_cut;
	JSINFO << cmd;
	VarMap var_map{};
  	po::store(po::command_line_parser(tokenize(cmd))
      .options(all_opts).positional(positional_opts).run(), var_map);
	TrentoGen_ = std::make_shared<trento::Collider>(var_map);
   	SetRanges(xymax, xymax, etamax);
	SetSteps(dxy, dxy, deta);
	JSINFO << "TRENTo set";
}

bool compare_E(trento::records r1, trento::records r2) { return r1.mult > r2.mult; }

std::pair<double, double> TrentoInitial::GenCenTab(std::string proj, std::string targ, VarMap var_map, int cL, int cH) {
	// Terminate for nonsense
	if (cL<0 || cL >100 || cH<0 || cH >100 || cH < cL) {
		JSWARN << "Wrong centrality cuts! To be terminated.";
		exit(-1);
	}
	// These are all the parameters that could change the shape of centrality tables
	// Normalization prefactor parameter is factorized
	// They form a table header
	trento::Collider another_collider(var_map);
	double beamE = var_map["beam-energy"].as<double>();
	double xsection = var_map["cross-section"].as<double>();
	double pvalue = var_map["reduced-thickness"].as<double>();
	double fluct = var_map["fluctuation"].as<double>();
	double nuclw = var_map["nucleon-width"].as<double>();
	double dmin = var_map["nucleon-min-dist"].as<double>();
	char buffer[512];
	std::sprintf(buffer, "%s-%s-E-%1.0f-X-%1.2f-p-%1.2f-k-%1.2f-w-%1.2f-d-%1.2f", 
				proj.c_str(), targ.c_str(), beamE, xsection, pvalue, fluct, nuclw, dmin);
	std::string header(buffer);
	JSINFO << "TRENTO centrality table header: " << header;
    // Create headering string hash tage for these parameter combination
	// Use this tag as a unique table filename for this specific parameter set
	std::hash<std::string> hash_function;
	size_t header_hash = hash_function(header);
	JSINFO << "Hash tag for this header: " << header_hash;
	// create dir incase it does not exist
	std::system("mkdir -p ./trento_data");
	char filename[512];
	std::sprintf(filename, "./trento_data/%zu", header_hash);
	// Step1: check it a table exist
	std::ifstream infile(filename);
	double Etab[101];
	double buff1, buff2;
	std::string line;
	if (infile.good()) {
		JSINFO << "The required centrality table exists. Load the table."; int i=0;
		while (std::getline(infile, line)) {
			if(line[0] != '#'){
				std::istringstream iss(line);
				iss >> buff1 >> buff2 >> Etab[i];
				i++;
			}
		}
		infile.close();
	}
	else {
		JSINFO << "TRENTo is generating new centrality table for this new parameter set";
		JSINFO << "It may take 10(s) to 1(min).";

		another_collider.run_events();
		// Get all records and sort according to totoal energy
		auto event_records = another_collider.all_records();
		std::sort(event_records.begin(), event_records.end(), compare_E);
		// write centrality table
		int nstep = std::ceil(event_records.size()/100);
		std::ofstream fout(filename);
		fout << "#\tproj\ttarj\tsqrts\tx\tp\tk\tw\td\n"
			 << "#\t" << proj << "\t" << targ << "\t"
			 << beamE  << "\t"
			 << xsection  << "\t"
			 << pvalue << "\t"
			 << fluct << "\t"
			 << nuclw << "\t"
			 << dmin << "\n"
			 << "#\tcen_L\tcen_H\tun-normalized total density\n";
		Etab[0] = 1e10;
		for (int i=1; i<100; i+=1) {
			auto ee = event_records[i*nstep];
			fout << i-1 << "\t" << i << "\t" << ee.mult << std::endl;
			Etab[i] = ee.mult;
		}
		auto ee = event_records.back();
		fout << 99 << "\t" << 100 << "\t" << ee.mult << std::endl;
		Etab[100] = ee.mult;
		fout.close();
	}
	JSINFO << "#########" << Etab[cL] << " " << Etab[cH];
	return std::make_pair(Etab[cL], Etab[cH]);
}

void TrentoInitial::Exec() {  
	JSINFO << " Exec TRENTo initial condition ";
	TrentoGen_->run_events();

	JSINFO << " TRENTo event info: ";
	auto tmp_event = TrentoGen_->expose_event();
    	info_.impact_parameter = TrentoGen_->all_records().back().b;
    	info_.num_participant = tmp_event.npart();
	info_.num_binary_collisions = tmp_event.ncoll();
    	info_.total_entropy = tmp_event.multiplicity();
   	info_.ecc = tmp_event.eccentricity();
	info_.psi = tmp_event.participant_plane();
    	info_.xmid = -GetXMax()+tmp_event.mass_center_index().first*tmp_event.dxy();
	info_.ymid = -GetYMax()+tmp_event.mass_center_index().second*tmp_event.dxy();
	JSINFO << "b\tnpart\tncoll\tET\t(x-com, y-com) (fm)";
	JSINFO << info_.impact_parameter << "\t" 
		 << info_.num_participant << "\t" 
		 << info_.num_binary_collisions << "\t"
		 << info_.total_entropy << "\t"
		 << "("<< info_.xmid << ", " << info_.ymid << ")";

    	JSINFO << " Load TRENTo density and ncoll density to JETSCAPE memory ";
	auto density_field = tmp_event.density_grid();
	auto ncoll_field = tmp_event.TAB_grid();
	JSINFO << density_field.num_elements() << " density elements";
    	for (int i=0; i<density_field.num_elements(); i++) {
       		entropy_density_distribution_.push_back(density_field.data()[i]);
    	}
	JSINFO << ncoll_field.num_elements() << " ncoll elements";
    	for (int i=0; i<ncoll_field.num_elements(); i++) {
       		num_of_binary_collisions_.push_back(ncoll_field.data()[i]);
    	}
	JSINFO << " TRENTO event generated and loaded ";
}

void TrentoInitial::Clear() {
    VERBOSE(2) << " : Finish creating initial condition ";
    entropy_density_distribution_.clear();
    num_of_binary_collisions_.clear();
}



} // end namespace Jetscape
