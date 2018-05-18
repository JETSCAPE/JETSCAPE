#ifndef Rate_H
#define Rate_H

#include <cstdlib>
#include <vector>
#include <string>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include "StochasticBase.h"
#include "Xsection.h"

// N1: dimension of rate table, N2: dimension of Xsection table, 
// F: matrix element function type
template <size_t N1, size_t N2, typename F>
class Rate: public virtual StochasticBase<N1>{
private:
	std::shared_ptr<Xsection<N2, F>> X;
    scalar find_max(std::vector<double> parameters);
	scalar calculate_scalar(std::vector<double> parameters);
	fourvec calculate_fourvec(std::vector<double> parameters);
	tensor calculate_tensor(std::vector<double> parameters);
	double _mass, _degen;
	bool _active;
public:
	Rate(std::string Name, std::string configfile, F f);
	void sample(std::vector<double> arg, 
				std::vector< fourvec > & IS);
	void initX(std::string fname){X->init(fname);}
	void loadX(std::string fname){X->load(fname);}
	bool IsActive(void) {return _active;}
};

#endif
