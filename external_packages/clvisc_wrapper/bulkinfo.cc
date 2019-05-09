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

#include "include/Config.h"
#include "include/opencl_backend.h"
#include "include/bulkinfo.h"

namespace clvisc {


BulkInfo::BulkInfo(int nx, int ny, int neta,
           int nx_skip, int ny_skip, int neta_skip,
           const OpenclBackend & backend, 
           const std::string & compile_option): nx_(nx), ny_(ny), neta_(neta), 
           nx_skip_(nx_skip), ny_skip_(ny_skip), neta_skip_(neta_skip),
           backend_(backend), compile_option_(compile_option){
    nx_out_ = int(floor((nx_-1)/nx_skip_)) + 1;
    ny_out_ = int(floor((ny_-1)/ny_skip_)) + 1;
    neta_out_ = int(floor((neta_-1)/neta_skip_)) + 1;
    /** tell jetscape what data is stored 
    * in std::vector<float> bulk_data_ */
    data_info_.push_back("energy_density");
    data_info_.push_back("entropy_density");
    data_info_.push_back("temperature");
    data_info_.push_back("pressure");
    data_info_.push_back("vx");
    data_info_.push_back("vy");
    data_info_.push_back("vz");
    data_info_.push_back("qgp_fraction");
    data_info_.push_back("pi00");
    data_info_.push_back("pi01");
    data_info_.push_back("pi02");
    data_info_.push_back("pi03");
    data_info_.push_back("pi11");
    data_info_.push_back("pi12");
    data_info_.push_back("pi13");
    data_info_.push_back("pi22");
    data_info_.push_back("pi23");
    data_info_.push_back("pi33");

    size_t size_out = nx_out_ * ny_out_ * neta_out_ * data_info_.size();
    h_bulk3d_1step_.resize(size_out);
    try {
        // build kernels for hydrodynamic evolution
        auto prg = backend_.BuildProgram("clvisc_kernel/kernel_bulk3d.cl", compile_option_);
        kernel_bulk3d_ = cl::Kernel(prg, "kernel_bulk3d");
        d_bulk3d_1step_ = backend_.CreateBuffer(size_out * sizeof(float));
    } catch (cl::Error & err ){
        std::cerr<<"Error:"<<err.what()<<"("<<err.err()<<")\n";
        std::cerr<<"@" << __FILE__ << ":line " << __LINE__ << std::endl;
        std::cerr<<ErrorMessage(err.err())<<std::endl;
        throw(err);
    }
}

/* push the fluid cell info of one time step back into a big vector */
void BulkInfo::add_data(const cl::Buffer & d_ev, 
                        const cl::Buffer & d_shear_pi,
                        const cl::Buffer & d_bulk_pi,
                        const cl::Image2D & eos_table){
    kernel_bulk3d_.setArg(0, d_bulk3d_1step_);
    kernel_bulk3d_.setArg(1, d_ev);
    kernel_bulk3d_.setArg(2, d_shear_pi);
    kernel_bulk3d_.setArg(3, d_bulk_pi);
    kernel_bulk3d_.setArg(4, eos_table);
    backend_.enqueue_run(kernel_bulk3d_,
            cl::NDRange(nx_out_, ny_out_, neta_out_),
            cl::NullRange);
    backend_.enqueue_copy(d_bulk3d_1step_, h_bulk3d_1step_);
    for (auto val : h_bulk3d_1step_) {
        bulk_data_.push_back(val);
    }
}

const std::vector<float> & BulkInfo::get_data() {
    return bulk_data_;
}

const std::vector<std::string> & BulkInfo::get_data_info() {
    return data_info_;
}

void BulkInfo::save(const std::string & fpath){
  // save the data to given fpath
    std::ofstream fout(fpath);
  int idx = 0;
  fout << "# nx=" << nx_out_ << std::endl;
  fout << "# ny=" << ny_out_ << std::endl;
  fout << "# nz=" << neta_out_ << std::endl;
  fout << "# num_entries =" << data_info_.size() << std::endl;
  fout << "# ";
  for (auto val : data_info_) {
      fout << val << " ";
  }
  fout << std::endl;
  for (auto val : bulk_data_) {
      idx ++;
      fout << val << " ";
      if (idx % data_info_.size() == 0) {
          fout << std::endl;
      }
  }
  fout.close();
}

} // end namespace clvisc


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


