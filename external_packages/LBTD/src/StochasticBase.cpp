#include "StochasticBase.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <thread>
#include "simpleLogger.h"
template<size_t N>
StochasticBase<N>::StochasticBase(std::string Name, std::string configfile):
_Name(Name)
{
	// read configfile
	boost::property_tree::ptree config;
	std::ifstream input(configfile);
	read_xml(input, config);
	
	std::vector<std::string> strs, slots;
	boost::split(strs, _Name, boost::is_any_of("/"));
	auto model_name = strs[0];
	auto process_name = strs[1];
	auto quantity_name = strs[2];

	auto tree = config.get_child(model_name+"."+process_name+"."+quantity_name);
	std::string allslots = tree.get<std::string>("<xmlattr>.slots");
	boost::split(slots, allslots, boost::is_any_of(",") );

	std::vector<size_t> shape;
	std::vector<double> low, high;
	for(auto & v : slots){
		shape.push_back(tree.get<size_t>("N"+v));
		low.push_back(tree.get<double>("L"+v));
		high.push_back(tree.get<double>("H"+v));
	}
    _FunctionMax = 
		std::make_shared<TableBase<scalar, N>>(Name+"/fmax", shape, low, high);
	_ZeroMoment = 
		std::make_shared<TableBase<scalar, N>>(Name+"/scalar", shape, low, high);
	_FirstMoment = 
		std::make_shared<TableBase<fourvec, N>>(Name+"/vector", shape, low, high);
	_SecondMoment = 
		std::make_shared<TableBase<tensor, N>>(Name+"/tensor", shape, low, high);
}

template<size_t N>
void StochasticBase<N>::load(std::string fname){
	LOG_INFO << "Loading " << _Name+"/fmax";
    _FunctionMax->Load(fname);
	LOG_INFO << "Loading " << _Name+"/scalar";
	_ZeroMoment->Load(fname);
	LOG_INFO << "Loading " << _Name+"/vector";
	_FirstMoment->Load(fname);
	LOG_INFO << "Loading " << _Name+"/tensor";
	_SecondMoment->Load(fname);
}


template<size_t N>
void StochasticBase<N>::init(std::string fname){
	LOG_INFO << _Name << " Generating tables";
	auto code = [this](int start, int end) { this->compute(start, end); };
	std::vector<std::thread> threads;
	size_t nthreads = std::thread::hardware_concurrency();
	size_t padding = size_t(std::ceil(_ZeroMoment->length()*1./nthreads));
	for(auto i=0; i<nthreads; ++i) {
		int start = i*padding;
		int end = std::min(padding*(i+1), _ZeroMoment->length());
		threads.push_back( std::thread(code, start, end) );
	}
	for(auto& t : threads) t.join();
	_FunctionMax->Save(fname);
	_ZeroMoment->Save(fname);
	_FirstMoment->Save(fname);
	_SecondMoment->Save(fname);
}

template<size_t N>
void StochasticBase<N>::compute(int start, int end){
	std::vector<size_t> index;
	index.resize(N);
	for(auto i=start; i<end; ++i){
		size_t q = i;
		for(int d=N-1; d>=0; d--){
			size_t dim = _ZeroMoment->shape(d);	
			size_t n = q%dim;
			q = q/dim;
			index[d] = n;
		}
		_FunctionMax->SetTableValue(index, 
						find_max(_FunctionMax->parameters(index))	);
		_ZeroMoment->SetTableValue(index, 
						calculate_scalar(_ZeroMoment->parameters(index))	);
		_FirstMoment->SetTableValue(index, 
						calculate_fourvec(_FirstMoment->parameters(index))	);
		_SecondMoment->SetTableValue(index, 
						calculate_tensor(_SecondMoment->parameters(index))	);
	}
}

template class StochasticBase<2>;
template class StochasticBase<3>;
template class StochasticBase<4>;
