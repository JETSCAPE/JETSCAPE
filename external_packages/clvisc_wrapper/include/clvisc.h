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


#ifndef __CL_VISC__
#define __CL_VISC__
/*!< Use cpp exceptionn to handel errors */
#define __CL_ENABLE_EXCEPTIONS 
// System includes
//#include <CL/cl.hpp>
#include "clideal.h"
#include "bulkinfo.h"

namespace clvisc {

/*! \class CLIdeal c++ wrapper for JetScape 
 *  */
class CLVisc
{
    private:
    Config cfg_;
    // CLVisc uses many data and member functions in ideal_
    CLIdeal ideal_;
    OpenclBackend backend_;
    double tau_;
    // size_ = cfg_.nx * cfg_.ny * cfg_.nz
    size_t size_;  
    std::string compile_option_;

    // notice h_share_pi_ has 10 components
    // pi00, 01, 02, 03, 11, 12, 13, 22, 23, 33
    std::vector<cl_real> h_shear_pi_;
    // notice h_bulk_pi_ has 1 component
    std::vector<cl_real> h_bulk_pi_;
    // notice h_net_charge_ has 4 components
    // net_baryon, net_electric, net_strangeness and one slot for future
    std::vector<cl_real4> h_net_charge_;

    cl::Buffer d_shear_pi_[3];
    cl::Buffer d_bulk_pi_[3];
    cl::Buffer d_net_charge_[3];

    // source term for IS equations
    cl::Buffer d_is_shear_pi_src_;
    cl::Buffer d_is_bulk_pi_src_;
    cl::Buffer d_udx_;
    cl::Buffer d_udy_;
    cl::Buffer d_udz_;
    cl::Buffer d_udiff_;

    // gpu buffer to check the goodness of each cell
    cl::Buffer d_goodcell_;

    // kernel functions in kernel_IS.cl
	cl::Kernel kernel_is_initialize_;
	cl::Kernel kernel_is_src_christoffel_;
	cl::Kernel kernel_is_src_alongx_;
	cl::Kernel kernel_is_src_alongy_;
	cl::Kernel kernel_is_src_alongz_;
	cl::Kernel kernel_is_update_pimn_;
	cl::Kernel kernel_is_get_udiff_;

    // kernel functions in kernel_visc.cl
	cl::Kernel kernel_visc_src_christoffel_;
	cl::Kernel kernel_visc_kt_src_alongx_;
	cl::Kernel kernel_visc_kt_src_alongy_;
	cl::Kernel kernel_visc_kt_src_alongz_;
	cl::Kernel kernel_visc_update_ev_;

    void initialize_gpu_buffer_();

    void israel_stewart_initialize_();

    // d_udiff = u_{visc}^{n} - u_{visc}^{n-1}
    void update_udiff_();

    void half_step_israel_stewart_(int step);
    void half_step_visc_(int step);
    // update half step using Runge-Kutta method, step = {1, 2}
    void half_step_(int step);

    public:

	CLVisc(const Config & cfg, std::string device_type, int device_id);

    BulkInfo bulkinfo_;

    inline const Config & get_config() {return cfg_;}

    // read initial energy density from external vector
    template <typename ValueType>
    void read_ini(const std::vector<ValueType> & ed);

    // read initial ed, vx, vy, vz vector
    template <typename ValueType>
    void read_ini(const std::vector<ValueType> & ed, 
                  const std::vector<ValueType> & vx, 
                  const std::vector<ValueType> & vy, 
                  const std::vector<ValueType> & vz);

    // read initial ed, vx, vy, vz vector and shear viscosity
    template <typename ValueType>
    void read_ini(const std::vector<ValueType> & ed, 
                  const std::vector<ValueType> & vx, 
                  const std::vector<ValueType> & vy, 
                  const std::vector<ValueType> & vz,
                  const std::vector<ValueType> & pi00,
                  const std::vector<ValueType> & pi01,
                  const std::vector<ValueType> & pi02,
                  const std::vector<ValueType> & pi03,
                  const std::vector<ValueType> & pi11,
                  const std::vector<ValueType> & pi12,
                  const std::vector<ValueType> & pi13,
                  const std::vector<ValueType> & pi22,
                  const std::vector<ValueType> & pi23,
                  const std::vector<ValueType> & pi33);

    // read initial ed, vx, vy, vz vector and shear viscosity,
    // bulk viscosity and charge current
    template <typename ValueType>
    void read_ini(const std::vector<ValueType> & ed, 
                  const std::vector<ValueType> & vx, 
                  const std::vector<ValueType> & vy, 
                  const std::vector<ValueType> & vz,
                  const std::vector<ValueType> & pi00,
                  const std::vector<ValueType> & pi01,
                  const std::vector<ValueType> & pi02,
                  const std::vector<ValueType> & pi03,
                  const std::vector<ValueType> & pi11,
                  const std::vector<ValueType> & pi12,
                  const std::vector<ValueType> & pi13,
                  const std::vector<ValueType> & pi22,
                  const std::vector<ValueType> & pi23,
                  const std::vector<ValueType> & pi33,
                  const std::vector<ValueType> & bulk_pi,
                  const std::vector<ValueType> & net_charge_baryon_,
                  const std::vector<ValueType> & net_charge_electric,
                  const std::vector<ValueType> & net_charge_strange);


    // run clvisc evolution for one time step
    void one_step();

    // run clvisc evolution for all time steps
    // stop when max_T < freeze_out_temperature
    void evolve();


	~CLVisc();

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


