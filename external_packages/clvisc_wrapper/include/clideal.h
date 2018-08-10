/*******************************************************************************
 * Copyright (c) 2018-2019 LongGang Pang, lgpang@qq.com
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and/or associated documentation files (the
 * "Materials"), to deal in the Materials without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Materials, and to
 * permit persons to whom the Materials are furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Materials.
 *
 * THE MATERIALS ARE PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * MATERIALS OR THE USE OR OTHER DEALINGS IN THE MATERIALS.
 ******************************************************************************/


#ifndef __CL_IDEAL__
#define __CL_IDEAL__
/*!< Use cpp exceptionn to handel errors */
#define __CL_ENABLE_EXCEPTIONS 
// System includes
//#include <CL/cl.hpp>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <ctime>
#include <algorithm>
#include <map>

#include <random>

#include "Config.h"
#include "opencl_backend.h"

namespace clvisc {

typedef struct
{
    int block_size; // num of threads along one dim on gpu
    int nx;  // number of grids along x
    int ny;  // number of grids along y
    int nz;  // number of grids along etas
    double tau0; // initial thermalization time
    double dt;   // time step
    float dx;   // x step 
    float dy;   // y step 
    float dz;   // etas step
    float etaos_xmin; // parameterized eta/s (T for minimum etaos)
    float etaos_ymin; // parameterized eta/s (minumum etaos)
    float etaos_left_slop; // parameterized eta/s (left slop)
    float etaos_right_slop; // parameterized eta/s (left slop)
    std::string result_directory;
} Config;

/*! \class CLIdeal c++ wrapper for JetScape 
 *  */
class CLIdeal
{
    private:
    double tau_;
    Config cfg_;
    OpenclBackend backend_;

    std::string compile_option_;
    // d_submax is used to compute the maximum
    // energy density of the fluctuating QGP
    cl::Buffer d_submax_;

    // stores the maximum energy density history
    std::vector<cl_real> max_ed_history_;

	cl::Kernel kernel_kt_src_christoffel_;
	cl::Kernel kernel_kt_src_alongx_;
	cl::Kernel kernel_kt_src_alongy_;
	cl::Kernel kernel_kt_src_alongz_;
	cl::Kernel kernel_update_ev_;
	cl::Kernel kernel_reduction_;

    void read_eos_table_(std::string fname, CompileOption &opts_);

    void initialize_gpu_buffer_();

     // update half step using Runge-Kutta method, step = {1, 2}
    void half_step_(int step);

    public:
    // h_ev_, d_ev_, eos_table_ will be used in class CLVisc
    // it would be good to make them public
    std::vector<cl_real4> h_ev_;
    cl::Buffer d_ev_[3];
    cl::Buffer d_src_;
    // image2d_t for eos_table
    cl::Image2D eos_table_;

	CLIdeal(const Config & cfg, std::string device_type, int device_id);

    // read initial energy density from external vector
    template <typename ValueType>
    void read_ini(const std::vector<ValueType> & ed);

    // read initial ed, vx, vy, vz vector
    template <typename ValueType>
    void read_ini(const std::vector<ValueType> & ed, 
                  const std::vector<ValueType> & vx, 
                  const std::vector<ValueType> & vy, 
                  const std::vector<ValueType> & vz);

    // run hydrodynamic evolution for one time step
    void one_step();

    // predict the first step to get u^{mu} for viscous hydro
    void predict_first_step();

    // return the maximum energy density, 
    float max_energy_density();

    // run hydrodynamic evolution for all time steps
    // stop when max_T < freeze_out_temperature
    void evolve();

    OpenclBackend & get_backend();
    std::string & get_compile_option();

	~CLIdeal();

};

} // end namespace clvisc

#endif 


/*! \mainpage 
 *  \section  */

/*! 
 *  \example 
*
 *  \code
 *
 *


 * \endcode
*/


