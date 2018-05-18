#ifndef MINIMIZER_H
#define MINIMIZER_H
#include <iostream>
#include <cmath>

#include <functional>
#include <memory>
#include <utility>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include "sampler.h"
#include "simpleLogger.h"

template < typename F >
double minimize_1d(F f, std::pair<double,double> const& range,
					double eps=0.01, int max_iter = 10, int N=5){
	double xlow = range.first, xhigh = range.second;
  	int loc=0, iter=0;
	double interval = xhigh-xlow;
	double value = 0.;
	double x, y, dx;
	do{
		iter ++;
		value = f(xlow);
		dx = (xhigh-xlow)/(N-1);
		for(int i=0; i< N; ++i){
			x = xlow+i*dx;
			y = f(x);
			if (value > y) {
				value = y;
				loc = i;
			}
		}
		xlow = std::max(xlow+dx*(loc-1), xlow);
		xhigh = std::min(xlow+dx*(loc+1), xhigh);
	}while(std::abs(xlow-xhigh)/interval > eps && iter < max_iter);
	value = f((xlow+xhigh)/2.);
    return value;
}

// GSL minimize N-dimensional
template <typename F>
class gsl_minimize_nd{
  F f;
  size_t dim;
  gsl_vector* ss;// starting x and step size 
  gsl_vector* x;
  std::unique_ptr <gsl_multimin_fminimizer,
					std::function<void(gsl_multimin_fminimizer*)>
					> workspace;
  

 
  static double f_wrapper(const gsl_vector * v, void * p){
	gsl_minimize_nd * t = reinterpret_cast<gsl_minimize_nd*>(p);
	double * x = new double[t->dim];
	for (int i=0; i< t->dim; ++i){
		x[i] = gsl_vector_get(v, i);
	}
	double val = t->f(x);
	delete [] x;
	return val;
  }

public:
  gsl_minimize_nd(F f, const size_t dim): 
    f(f), dim(dim),
    ss(gsl_vector_alloc(dim)),
    x(gsl_vector_alloc(dim)),
    workspace(
		gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, dim),
		gsl_multimin_fminimizer_free	)
	{}
  ~gsl_minimize_nd(){
	gsl_vector_free(ss);
	gsl_vector_free(x);
  }

  double minimize(std::vector<double> start, std::vector<double> step,
				  int max_iter, double eps){
	gsl_multimin_function minex_func;
	size_t iter = 0;
  	int status;
  	double size; // simplex size
	for (int i=0; i<dim; ++i){
		gsl_vector_set(ss, i, step[i]);
		gsl_vector_set(x, i, start[i]);
	}
	minex_func.n = dim;
	minex_func.f = f_wrapper;
	minex_func.params = this;
	gsl_multimin_fminimizer_set(workspace.get(), &minex_func, x, ss);
	do{
      iter++;
      status = gsl_multimin_fminimizer_iterate(workspace.get());
      if (status) break;
      size = gsl_multimin_fminimizer_size(workspace.get());
      status = gsl_multimin_test_size(size, eps);
    }while (status == GSL_CONTINUE && iter < max_iter);
    
    return (workspace.get())->fval;
  }
};

template <typename F>
double minimize_nd(F func, const size_t dim, 
					std::vector<double> start, std::vector<double> step,
                    int max_iter=100, double eps=0.01){
  return gsl_minimize_nd<F>(func, dim).minimize(start, step, max_iter, eps);
}

// Use MC random walk as a maximizer
// ----------Affine-invariant metropolis sample-------------------

template<typename F>
std::vector<double> MC_maximize(F f_, int n_dims_,
		std::vector<std::pair<double,double>> const& range, int steps=20){
	auto a = AiMS<F>(f_, n_dims_);
	a.sample(range, steps);
	return a.getmaxloc();
}

#endif 
