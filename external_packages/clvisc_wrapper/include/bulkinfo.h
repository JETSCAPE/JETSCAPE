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


#ifndef __CL_BULKINFO__
#define __CL_BULKINFO__
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

/*! \class Output module to save bulkinfo
 *  */
class BulkInfo
{
    private:
    int nx_, ny_, neta_;
    int nx_skip_;
    int ny_skip_;
    int neta_skip_;
    OpenclBackend backend_;
    std::string compile_option_;

	cl::Kernel kernel_bulk3d_;
    // store the medium information of whole evolution
    std::vector<float> bulk_data_;
    std::vector<std::string> data_info_;
    // collect 3-dimensional bulk info for 1 step
    cl::Buffer d_bulk3d_1step_;
    std::vector<float> h_bulk3d_1step_;
    // num of cells to output along each dimension
    int nx_out_;
    int ny_out_;
    int neta_out_;

    public:
    // using sparse lattice in the output to save memory
	BulkInfo(int nx, int ny, int nz,
           int nx_skip, int ny_skip, int neta_skip,
           const OpenclBackend & backend, 
           const std::string & compile_option);

    void add_data(const cl::Buffer & d_ev, 
                  const cl::Buffer & d_shear_pi,
                  const cl::Buffer & d_bulk_pi,
                  const cl::Image2D & eos_table);

    void save(const std::string & fpath);

    const std::vector<float> & get_data();

    const std::vector<std::string> & get_data_info();
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


