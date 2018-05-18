#ifndef Xsection_H
#define Xsection_H

#include <cstdlib>
#include <vector>
#include <string>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include "fast_exp.h"
#include "StochasticBase.h"


template <size_t N, typename F>
class Xsection: public virtual StochasticBase<N> {
private:
    scalar find_max(std::vector<double> parameters);
	scalar calculate_scalar(std::vector<double> parameters);
	fourvec calculate_fourvec(std::vector<double> parameters);
	tensor calculate_tensor(std::vector<double> parameters);
	double _mass;
	F _f;// the matrix element
	const FastExp<double> fast_exp_;
public:
	Xsection(std::string Name, std::string configfile, F f);
	void sample(std::vector<double> arg, 
						std::vector< fourvec > & FS);
};

#endif
