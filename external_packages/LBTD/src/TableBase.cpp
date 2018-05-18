#include "TableBase.h"
#include "H5Cpp.h"
#include "predefine.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include "simpleLogger.h"

// Default approximation function

template <typename T>
T default_approximate_function(Dvec values){
	return T::unity();
}

template <typename T, size_t N>
TableBase<T, N>::TableBase(std::string Name, Svec shape, Dvec low, Dvec high):
_Name(Name), _rank(N), _power_rank(std::pow(2, _rank)), 
_shape(shape), _low(low), _high(high),_table(_shape)
{
	LOG_INFO<<_Name << " dim=" << _rank;
	for(auto i=0; i<_rank; ++i){
		_step.push_back((high[i]-low[i])/(shape[i]-1));
	}
	// Set default approximation function to return 1
	ApproximateFunction = default_approximate_function<T>;
}

template <typename T, size_t N>
T TableBase<T, N>::InterpolateTable(Dvec values){
   Svec start_index;
   Dvec w;
   for(auto i=0; i<_rank; ++i) {
       auto x = (values[i]-_low[i])/_step[i];
       x = std::min(std::max(x, 0.), _shape[i]-2.); // cut at lower and higher bounds bounds
       size_t nx = size_t(std::floor(x));
       double rx = x-nx;
       start_index.push_back(nx);
       w.push_back(rx);
   }
   Svec index(_rank);
   T result{0.};
   Dvec corner_values(_rank); // hold x values at the corner of the hyper cube
   for(auto i=0; i<_power_rank; ++i) {
        auto W = 1.0;
        for (auto j=0; j<_rank; ++j) {
            index[j] = start_index[j] + ((i & ( 1 << j )) >> j);
            W *= (index[j]==start_index[j])?(1.-w[j]):w[j];
            corner_values[j] = _low[j] + _step[j]*index[j];
        }
        // We interp f/f_approx
        result = result + _table(index)/ApproximateFunction(corner_values)*W;
   }
   // multiply the interp function back with f_approx
   return result*ApproximateFunction(values);
}

template <typename T, size_t N>
void TableBase<T, N>::SetTableValue(Svec index, T v){
    _table(index) = v;
}

template <typename T, size_t N>
bool TableBase<T, N>::Save(std::string fname){
	H5::Exception::dontPrint(); // suppress error messages
	// the is a dumb implementation
	H5::H5File file;
LOG_INFO << "save start";
	if( boost::filesystem::exists(fname)) file = H5::H5File(fname, H5F_ACC_RDWR);
	else file = H5::H5File(fname, H5F_ACC_TRUNC);
LOG_INFO << "file opened";
	std::vector<std::string> levels;
	boost::split(levels, _Name, boost::is_any_of("/"));
	std::string prefix = "";
	H5::Group group;
	for (auto& v : levels) {
		prefix += ("/"+v);
  		try{
    		group = file.openGroup(prefix.c_str());
			if (prefix == "/"+_Name) {	// if the last group existed before
				// It need to be deleted and rebuild
				H5Ldelete(file.getId(), prefix.c_str(), H5P_DEFAULT);
				LOG_WARNING<<"old data deleted and will be overwirtten";
				group = file.createGroup(prefix.c_str());
			}
  		}catch (...) {
    		group = file.createGroup(prefix.c_str());
 		}
	}
LOG_INFO << "Create group";
	hdf5_add_scalar_attr(group, "rank", _rank);
	for (auto i=0; i<_rank; ++i){
		hdf5_add_scalar_attr(group, "shape-"+std::to_string(i), _shape[i]);
		hdf5_add_scalar_attr(group, "low-"+std::to_string(i), _low[i]);
		hdf5_add_scalar_attr(group, "high-"+std::to_string(i), _high[i]);
	}
	
	boost::multi_array<double, N> buffer(_shape);
	hsize_t dims[_rank];
	for (auto i=0; i<_rank; ++i) dims[i]=_shape[i];
	H5::DSetCreatPropList proplist{};
	proplist.setChunk(_rank, dims);

	H5::DataSpace dataspace(_rank, dims);
	auto datatype(H5::PredType::NATIVE_DOUBLE);
	
	for(auto comp=0; comp<T::size(); ++comp) {
		for(auto i=0; i<_table.num_elements(); ++i) {
			T item = _table.data()[i];
			buffer.data()[i] = item.get(comp);
		}
		auto dsname = prefix+"/"+std::to_string(comp);
		H5::DataSet dataset = file.createDataSet(dsname, datatype, dataspace, proplist);
		dataset.write(buffer.data(), datatype);
	}
	file.close();
	return true;
}

template <typename T, size_t N>
bool TableBase<T, N>::Load(std::string fname){
	H5::H5File file(fname, H5F_ACC_RDONLY);
	H5::Group group = H5::Group( file.openGroup( "/"+_Name ));
	size_t temp_rank;
	hdf5_read_scalar_attr(group, "rank", temp_rank);
	if (temp_rank != _rank) {
		LOG_FATAL<< "Table rank does not match";
		file.close();
		return false;
	}
	else{
		LOG_INFO<< "Rank compitable, loading table";
		for (auto i=0; i<_rank; ++i){
			hdf5_read_scalar_attr(group, "shape-"+std::to_string(i), _shape[i]);
			hdf5_read_scalar_attr(group, "low-"+std::to_string(i), _low[i]);
			hdf5_read_scalar_attr(group, "high-"+std::to_string(i), _high[i]);
			_step[i] = (_high[i] - _low[i])/(_shape[i]-1.);
		}
		_table.resize(_shape);
		hsize_t dims[_rank];
		for (auto i=0; i<_rank; ++i) dims[i]=_shape[i];
		boost::multi_array<double, N> buffer(_shape);
		H5::DataSpace dataspace(_rank, dims);
		auto datatype(H5::PredType::NATIVE_DOUBLE);
		for(auto comp=0; comp<T::size(); ++comp) {
			H5::DataSet dataset = file.openDataSet("/"+_Name+"/"+std::to_string(comp));
			dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE,
						 dataspace, dataset.getSpace());
			for(auto i=0; i<_table.num_elements(); ++i) {
				_table.data()[i].set(comp, buffer.data()[i]);
			}
		}
		file.close();
	}
	return true;
}

template class TableBase<scalar, 2>;
template class TableBase<scalar, 3>;
template class TableBase<scalar, 4>;
template class TableBase<fourvec, 2>;
template class TableBase<fourvec, 3>;
template class TableBase<fourvec, 4>;
template class TableBase<tensor, 2>;
template class TableBase<tensor, 3>;
template class TableBase<tensor, 4>;



