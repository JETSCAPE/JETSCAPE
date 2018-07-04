#ifndef TABLE_BASE_H
#define TABLE_BASE_H

#include <vector>
#include <string>
#include <boost/multi_array.hpp>
#include <iostream>
#include "lorentz.h"

typedef std::vector<double> Dvec;
typedef std::vector<size_t> Svec;

// Base class of a table of type T with dimension N
template <typename T, size_t N>
class TableBase{
protected:
    const std::string _Name;
    const size_t _rank, _power_rank;
    Svec _shape;
    Dvec _low, _high;
    Dvec _step;
    boost::multi_array<T, N> _table;
    T(*ApproximateFunction)(Dvec values);
public:
	TableBase(std::string, Svec, Dvec, Dvec);
	T InterpolateTable(Dvec values);
    void SetTableValue(Svec index, T v);
    void SetApproximateFunction(T(*f)(Dvec values)){
    	ApproximateFunction = f;
    	};
    bool Save(std::string);
    bool Load(std::string);
	size_t shape(size_t i) {return _shape[i];}
	size_t rank(void) {return _rank;}
	size_t length(void) {
		size_t result=1; 
		for(auto& D : _shape) result *= D;
		return result;
	}
	Dvec parameters(Svec index){
		Dvec res;
		for(size_t i=0; i<_rank; i++)
			res.push_back(_low[i] + _step[i]*index[i]);
		return res;
	}
};





#endif

