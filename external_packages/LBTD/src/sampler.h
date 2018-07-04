#ifndef SAMPLER_H
#define SAMPLER_H

#include <iostream>
#include <cmath>
#include <functional>
#include <memory>
#include <utility>
#include "random.h"
#include "simpleLogger.h"
#include "stat.h"

template < typename F >
double sample_1d(F f, std::pair<double,double> const& range, double fmax){
	int limit = 10000;
  	double y, x, xlow=range.first, xhigh=range.second;
  	double interval = xhigh-xlow;
  	int counter = 0;
	do{
		x = xlow+Srandom::init_dis(Srandom::gen)*interval;
		y = f(x)/fmax;
		if (y > 1.0) LOG_WARNING << "1d rejection, f/fmax = " << y << " > 1";
		counter ++;
	}while(Srandom::rejection(Srandom::gen)>y && counter < limit);
	if(counter==limit) LOG_WARNING <<  "1d rejection, too many tries = " << limit;
	SamplerStat::count_1d ++; SamplerStat::total_1d += counter;
	return x;
}

template < typename F >
std::vector<double> sample_nd(F f, int dim, std::vector<std::pair<double,double>> const& range, double fmax, bool & status){
	int limit = 50000;
	double * x = new double[dim];
  	double y;
  	double * interval = new double[dim];
	int counter = 0;
	for(int i=0; i<dim; i++) interval[i] = range[i].second - range[i].first;
	do{
		// random choice
		for(int i=0; i<dim; i++) 
			x[i] = range[i].first+Srandom::init_dis(Srandom::gen)*interval[i];
		y = f(x)/fmax;
		if (y > 1.0) {
			LOG_WARNING << "nd rejection, f/fmax = " << y << " > 1";
			status = false;
		}
		counter ++;
	}while(Srandom::rejection(Srandom::gen)>y && counter < limit);
	std::vector<double> res(dim);
	for(int i=0; i<dim; i++) res[i] = x[i];
	delete[] x;
	delete[] interval;
	if(counter==limit) {
		LOG_WARNING <<  "nd rejection, too many tries = " << limit;
		status = false;
	}
	SamplerStat::count_nd ++; SamplerStat::total_nd += counter;
	return res;
}

// ----------Affine-invariant metropolis sample-------------------
struct walker{
	double * posi;
	double P;
};

template < typename F >
class AiMS{
private:
	F f;
	size_t n_dims, Nwalker;
	std::vector<walker> walkers, buff_walkers;
	double maxP;
	std::vector<double> maxloc;
	
	void initialize(std::vector<std::pair<double,double>> range){
		for (size_t i=0; i<Nwalker; ++i){
			do{
				for (size_t j=0; j < n_dims; ++j){
					walkers[i].posi[j] = range[j].first 
							+ (range[j].second-range[j].first)
							  *Srandom::init_dis(Srandom::gen);
				}
				walkers[i].P = f(walkers[i].posi);
			} while(walkers[i].P <= 1e-22);
			for (size_t j=0; j < n_dims; ++j)
					buff_walkers[i].posi[j] = walkers[i].posi[j];
			buff_walkers[i].P = walkers[i].P;
		}
	}
	void update(void){
		size_t ri;
		double sqz, z, Ptry, Paccept;
		double * xtry = new double[n_dims];
		walker w, wr;
		for (size_t i=0; i<Nwalker; ++i){
			do{ 
				ri = std::rand() % Nwalker;
			}while(i==ri);
			w = walkers[i];
			wr = walkers[ri];
			sqz = Srandom::sqrtZ(Srandom::gen);
			z = sqz*sqz;
			for (size_t j=0; j < n_dims; ++j) 
				xtry[j] = wr.posi[j] + z*(w.posi[j] - wr.posi[j]);
			Ptry = f(xtry); 
			// A side product is to find maximum in a Monte Carlo way
			if (Ptry > maxP){
				maxP = Ptry;
				for (size_t j=0; j < n_dims; ++j) maxloc[j] = xtry[j];
			}
			Paccept = Ptry/w.P*std::pow(z, n_dims-1);
			if (Paccept >= 1.0){
				for (size_t j=0; j < n_dims; ++j) buff_walkers[i].posi[j] = xtry[j];
				buff_walkers[i].P = Ptry;
			}
			else if (Paccept >= Srandom::rejection(Srandom::gen)){
				for (size_t j=0; j < n_dims; ++j) buff_walkers[i].posi[j] = xtry[j];
				buff_walkers[i].P = Ptry;
			}
		}
		for (size_t i=0; i<Nwalker; ++i){
			for (size_t j=0; j < n_dims; ++j) 
				walkers[i].posi[j] = buff_walkers[i].posi[j];
			walkers[i].P = buff_walkers[i].P;
		}
		delete[] xtry;
	}

public:
	AiMS(F f_, int n_dims_): f(f_), n_dims(n_dims_), Nwalker(n_dims*4){
		maxloc.resize(n_dims);
		walkers.resize(Nwalker); 
		buff_walkers.resize(Nwalker);
		for (auto&& w : walkers) w.posi = new double[n_dims];
		for (auto&& w : buff_walkers) w.posi = new double[n_dims];
	}
	~AiMS(){
		for (auto&& w : walkers) delete[] w.posi;
		for (auto&& w : buff_walkers) delete[] w.posi;
	}
	std::vector<double> sample(std::vector<std::pair<double,double>> const& range_, int steps_){
		initialize(range_);
		maxP = 0.;
		for (size_t i = 0; i<steps_; i++) update();
		std::vector<double> result(n_dims);
		for (size_t i = 0; i<n_dims; i++) result[i] = walkers[0].posi[i];
		return result;
	}
	double getMax(void) {return maxP;}
	std::vector<double> getmaxloc(void) {return maxloc;}
};

template<typename F>
std::vector<double> MC_sample(F f_, int n_dims_,
		std::vector<std::pair<double,double>> const& range, int steps=20){
	return AiMS<F>(f_, n_dims_).sample(range, steps);
}
#endif 
