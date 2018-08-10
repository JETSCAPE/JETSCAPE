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

#include <cmath>
#include <cstring>
#include <ctime>
#include <algorithm>
#include "include/clvisc.h"

namespace clvisc {
CLVisc::CLVisc(const Config & cfg, std::string device_type,
        int device_id):cfg_(cfg),
        ideal_(CLIdeal(cfg, device_type, device_id)),
        backend_(ideal_.get_backend())
{
    // set current time to initial time
    tau_ = cfg.tau0;

    size_ = cfg.nx * cfg.ny * cfg.nz;

    compile_option_ = ideal_.get_compile_option();

    std::cout << compile_option_ << std::endl;

    try {
        // build kernels for hydrodynamic evolution
        auto prg = backend_.BuildProgram("clvisc_kernel/kernel_IS.cl",
                compile_option_);
        // is == Israel Stewart here
        kernel_is_initialize_ = cl::Kernel(prg, "visc_initialize");
        kernel_is_src_christoffel_ = cl::Kernel(prg, "visc_src_christoffel");
        kernel_is_src_alongx_ = cl::Kernel(prg, "visc_src_alongx");
        kernel_is_src_alongy_ = cl::Kernel(prg, "visc_src_alongy");
        kernel_is_src_alongz_ = cl::Kernel(prg, "visc_src_alongz");
        kernel_is_update_pimn_ = cl::Kernel(prg, "update_pimn");
        kernel_is_get_udiff_ = cl::Kernel(prg, "get_udiff");

        auto prg2 = backend_.BuildProgram("clvisc_kernel/kernel_visc.cl",
                compile_option_);
        kernel_visc_src_christoffel_ = cl::Kernel(prg2, "kt_src_christoffel");
        kernel_visc_kt_src_alongx_ = cl::Kernel(prg2, "kt_src_alongx");
        kernel_visc_kt_src_alongy_ = cl::Kernel(prg2, "kt_src_alongy");
        kernel_visc_kt_src_alongz_ = cl::Kernel(prg2, "kt_src_alongz");
        kernel_visc_update_ev_ = cl::Kernel(prg2, "update_ev");
    } catch (cl::Error & err ){
        std::cerr<<"Error:"<<err.what()<<"("<<err.err()<<")\n";
    }
}


void CLVisc::initialize_gpu_buffer_() {
    bool read_only_option = false;
    // ed_bytes is the size of energy density in bytes
    // d_ev_[i] has size (4*ed_bytes)
    // d_shear_pi_[i] has size (10*ed_bytes)
    size_t ed_bytes = sizeof(cl_real) * size_;

    for (int i = 0; i < 3; i++) {
        d_shear_pi_[i] = backend_.CreateBufferByCopyVector(h_shear_pi_, read_only_option);
        d_bulk_pi_[i] = backend_.CreateBufferByCopyVector(h_bulk_pi_, read_only_option);
        d_net_charge_[i] = backend_.CreateBufferByCopyVector(h_net_charge_, read_only_option);
    }

    // initialize the d_udz_ to vector of (real4)
    // in case the nz = 1 and one wants to skip z direction calculation
    std::vector<cl_real4> h_udz_(size_, (cl_real4){0.0f, 0.0f, 0.0f, 0.0f});

    d_is_shear_pi_src_ = backend_.CreateBuffer(10*ed_bytes);
    d_is_bulk_pi_src_ = backend_.CreateBuffer(ed_bytes);
    d_udx_ = backend_.CreateBuffer(4*ed_bytes);    // cl_real4 type
    d_udy_ = backend_.CreateBuffer(4*ed_bytes);
    d_udz_ = backend_.CreateBufferByCopyVector(h_udz_, read_only_option);
    d_udiff_ = backend_.CreateBuffer(4*ed_bytes);
    d_goodcell_ = backend_.CreateBuffer(ed_bytes);
}


void CLVisc::israel_stewart_initialize_() {
    kernel_is_initialize_.setArg(0, d_shear_pi_[1]);
    kernel_is_initialize_.setArg(1, d_goodcell_);
    kernel_is_initialize_.setArg(2, d_udiff_);
    kernel_is_initialize_.setArg(3, ideal_.d_ev_[1]);
    kernel_is_initialize_.setArg(4, static_cast<cl_real> (tau_));
    kernel_is_initialize_.setArg(5, ideal_.eos_table_);
    backend_.enqueue_run(kernel_is_initialize_, cl::NDRange(size_), cl::NullRange);
}


void CLVisc::update_udiff_() {
    kernel_is_get_udiff_.setArg(0, d_udiff_);
    kernel_is_get_udiff_.setArg(1, ideal_.d_ev_[0]);
    kernel_is_get_udiff_.setArg(2, ideal_.d_ev_[1]);
    backend_.enqueue_run(kernel_is_get_udiff_,
            cl::NDRange(size_), cl::NullRange);
}

void CLVisc::half_step_israel_stewart_(int step) {
    kernel_is_src_christoffel_.setArg(0, d_is_shear_pi_src_);
    kernel_is_src_christoffel_.setArg(1, d_shear_pi_[step]);
    kernel_is_src_christoffel_.setArg(2, ideal_.d_ev_[step]);
    kernel_is_src_christoffel_.setArg(3, static_cast<cl_real> (tau_));
    kernel_is_src_christoffel_.setArg(4, step);
    backend_.enqueue_run(kernel_is_src_christoffel_,
            cl::NDRange(size_), cl::NullRange);

    kernel_is_src_alongx_.setArg(0, d_is_shear_pi_src_);
    kernel_is_src_alongx_.setArg(1, d_udx_);
    kernel_is_src_alongx_.setArg(2, d_shear_pi_[step]);
    kernel_is_src_alongx_.setArg(3, ideal_.d_ev_[step]);
    kernel_is_src_alongx_.setArg(4, ideal_.eos_table_);
    kernel_is_src_alongx_.setArg(5, static_cast<cl_real> (tau_));
    backend_.enqueue_run(kernel_is_src_alongx_, 
            cl::NDRange(cfg_.block_size, cfg_.ny, cfg_.nz),
            cl::NDRange(cfg_.block_size, 1, 1));

    kernel_is_src_alongy_.setArg(0, d_is_shear_pi_src_);
    kernel_is_src_alongy_.setArg(1, d_udy_);
    kernel_is_src_alongy_.setArg(2, d_shear_pi_[step]);
    kernel_is_src_alongy_.setArg(3, ideal_.d_ev_[step]);
    kernel_is_src_alongy_.setArg(4, ideal_.eos_table_);
    kernel_is_src_alongy_.setArg(5, static_cast<cl_real> (tau_));
    backend_.enqueue_run(kernel_is_src_alongy_, 
            cl::NDRange(cfg_.nx, cfg_.block_size, cfg_.nz),
            cl::NDRange(1, cfg_.block_size, 1));

    if (cfg_.nz != 1) {
        kernel_is_src_alongz_.setArg(0, d_is_shear_pi_src_);
        kernel_is_src_alongz_.setArg(1, d_udz_);
        kernel_is_src_alongz_.setArg(2, d_shear_pi_[step]);
        kernel_is_src_alongz_.setArg(3, ideal_.d_ev_[step]);
        kernel_is_src_alongz_.setArg(4, ideal_.eos_table_);
        kernel_is_src_alongz_.setArg(5, static_cast<cl_real> (tau_));
        backend_.enqueue_run(kernel_is_src_alongz_, 
                cl::NDRange(cfg_.nx, cfg_.ny, cfg_.block_size),
                cl::NDRange(1, 1, cfg_.block_size));
    }

    kernel_is_update_pimn_.setArg(0, d_shear_pi_[3-step]);
    kernel_is_update_pimn_.setArg(1, d_goodcell_);
    kernel_is_update_pimn_.setArg(2, d_shear_pi_[1]);
    kernel_is_update_pimn_.setArg(3, d_shear_pi_[step]);
    kernel_is_update_pimn_.setArg(4, ideal_.d_ev_[1]);
    kernel_is_update_pimn_.setArg(5, ideal_.d_ev_[2]);
    kernel_is_update_pimn_.setArg(6, d_udiff_);
    kernel_is_update_pimn_.setArg(7, d_udx_);
    kernel_is_update_pimn_.setArg(8, d_udy_);
    kernel_is_update_pimn_.setArg(9, d_udz_);
    kernel_is_update_pimn_.setArg(10, d_is_shear_pi_src_);
    kernel_is_update_pimn_.setArg(11, ideal_.eos_table_);
    kernel_is_update_pimn_.setArg(12, static_cast<cl_real> (tau_));
    kernel_is_update_pimn_.setArg(13, step);
    backend_.enqueue_run(kernel_is_update_pimn_,
            cl::NDRange(size_), cl::NullRange);

}

void CLVisc::half_step_visc_(int step) {
    kernel_visc_src_christoffel_.setArg(0, ideal_.d_src_);
    kernel_visc_src_christoffel_.setArg(1, ideal_.d_ev_[step]);
    kernel_visc_src_christoffel_.setArg(2, d_shear_pi_[step]);
    kernel_visc_src_christoffel_.setArg(3, ideal_.eos_table_);
    kernel_visc_src_christoffel_.setArg(4, static_cast<cl_real> (tau_));
    kernel_visc_src_christoffel_.setArg(5, step);
    backend_.enqueue_run(kernel_visc_src_christoffel_,
            cl::NDRange(size_), cl::NullRange);

    kernel_visc_kt_src_alongx_.setArg(0, ideal_.d_src_);
    kernel_visc_kt_src_alongx_.setArg(1, ideal_.d_ev_[step]);
    kernel_visc_kt_src_alongx_.setArg(2, d_shear_pi_[step]);
    kernel_visc_kt_src_alongx_.setArg(3, ideal_.eos_table_);
    kernel_visc_kt_src_alongx_.setArg(4, static_cast<cl_real> (tau_));
    backend_.enqueue_run(kernel_visc_kt_src_alongx_,
            cl::NDRange(cfg_.block_size, cfg_.ny, cfg_.nz),
            cl::NDRange(cfg_.block_size, 1, 1));

    kernel_visc_kt_src_alongy_.setArg(0, ideal_.d_src_);
    kernel_visc_kt_src_alongy_.setArg(1, ideal_.d_ev_[step]);
    kernel_visc_kt_src_alongy_.setArg(2, d_shear_pi_[step]);
    kernel_visc_kt_src_alongy_.setArg(3, ideal_.eos_table_);
    kernel_visc_kt_src_alongy_.setArg(4, static_cast<cl_real> (tau_));
    backend_.enqueue_run(kernel_visc_kt_src_alongy_,
            cl::NDRange(cfg_.nx, cfg_.block_size, cfg_.nz),
            cl::NDRange(1, cfg_.block_size, 1));


    if (cfg_.nz != 1) {
        kernel_visc_kt_src_alongz_.setArg(0, ideal_.d_src_);
        kernel_visc_kt_src_alongz_.setArg(1, ideal_.d_ev_[step]);
        kernel_visc_kt_src_alongz_.setArg(2, d_shear_pi_[step]);
        kernel_visc_kt_src_alongz_.setArg(3, ideal_.eos_table_);
        kernel_visc_kt_src_alongz_.setArg(4, static_cast<cl_real> (tau_));
        backend_.enqueue_run(kernel_visc_kt_src_alongz_,
                cl::NDRange(cfg_.nx, cfg_.ny, cfg_.block_size),
                cl::NDRange(1, 1, cfg_.block_size));
    }

    kernel_visc_update_ev_.setArg(0, ideal_.d_ev_[3-step]);
    kernel_visc_update_ev_.setArg(1, ideal_.d_ev_[1]);
    kernel_visc_update_ev_.setArg(2, d_shear_pi_[0]);
    kernel_visc_update_ev_.setArg(3, d_shear_pi_[3-step]);
    kernel_visc_update_ev_.setArg(4, ideal_.d_src_);
    kernel_visc_update_ev_.setArg(5, ideal_.eos_table_);
    kernel_visc_update_ev_.setArg(6, static_cast<cl_real> (tau_));
    kernel_visc_update_ev_.setArg(7, step);
    backend_.enqueue_run(kernel_visc_update_ev_,
            cl::NDRange(size_), cl::NullRange);
} // end half_step_visc_


void CLVisc::half_step_(int step) {
    half_step_israel_stewart_(step);
    half_step_visc_(step);
}


template <typename ValueType>
void CLVisc::read_ini(const std::vector<ValueType> & ed) {
    ideal_.read_ini(ed);
    h_shear_pi_ = std::vector<cl_real> (10*size_, 0.0);
    h_bulk_pi_ = std::vector<cl_real> (size_, 0.0);
    cl_real4 zero4 = (cl_real4){0.0, 0.0, 0.0, 0.0};
    h_net_charge_ = std::vector<cl_real4> (size_, zero4);
    initialize_gpu_buffer_();
}

template <typename ValueType>
void CLVisc::read_ini(const std::vector<ValueType> & ed, 
                      const std::vector<ValueType> & vx, 
                      const std::vector<ValueType> & vy, 
                      const std::vector<ValueType> & vz)
{
    ideal_.read_ini(ed, vx, vy, vz);
    h_shear_pi_ = std::vector<cl_real> (10*size_, 0.0);
    h_bulk_pi_ = std::vector<cl_real> (size_, 0.0);
    cl_real4 zero4 = (cl_real4){0.0, 0.0, 0.0, 0.0};
    h_net_charge_ = std::vector<cl_real4> (size_, zero4);
    initialize_gpu_buffer_();
}

template <typename ValueType>
void CLVisc::read_ini(const std::vector<ValueType> & ed, 
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
              const std::vector<ValueType> & pi33) {
    ideal_.read_ini(ed, vx, vy, vz);
    h_shear_pi_.clear();
    for (size_t i = 0; i < size_; i++) {
        h_shear_pi_.push_back(static_cast<cl_real>(pi00.at(i)));
        h_shear_pi_.push_back(static_cast<cl_real>(pi01.at(i)));
        h_shear_pi_.push_back(static_cast<cl_real>(pi02.at(i)));
        h_shear_pi_.push_back(static_cast<cl_real>(pi03.at(i)));
        h_shear_pi_.push_back(static_cast<cl_real>(pi11.at(i)));
        h_shear_pi_.push_back(static_cast<cl_real>(pi12.at(i)));
        h_shear_pi_.push_back(static_cast<cl_real>(pi13.at(i)));
        h_shear_pi_.push_back(static_cast<cl_real>(pi22.at(i)));
        h_shear_pi_.push_back(static_cast<cl_real>(pi23.at(i)));
        h_shear_pi_.push_back(static_cast<cl_real>(pi33.at(i)));
    }
    h_bulk_pi_ = std::vector<cl_real> (size_, 0.0);
    cl_real4 zero4 = (cl_real4){0.0, 0.0, 0.0, 0.0};
    h_net_charge_ = std::vector<cl_real4> (size_, zero4);
    initialize_gpu_buffer_();
}




void CLVisc::evolve() {
    int max_loops = 5000;
    float total_exec_time = 0.0;
    std::time_t timer1, timer2;
    std::time(&timer1);
    try {
        israel_stewart_initialize_();
        ideal_.predict_first_step();
        for (int loop = 0; loop < max_loops; loop++) {
            std::cout << "tau = " << tau_ << " fm; " << std::endl;
            backend_.enqueue_copy(d_shear_pi_[1], d_shear_pi_[0], 10*size_*sizeof(cl_real));
            backend_.enqueue_copy(ideal_.d_ev_[1], ideal_.d_ev_[0], size_*sizeof(cl_real4));
            half_step_(1);
            tau_ += cfg_.dt;
            half_step_(2);
            update_udiff_();
            if (loop % 10 == 0) {
                float max_ed = ideal_.max_energy_density();
                if ( max_ed < 0.5 ) break;
                std::cout << "max_ed = " << max_ed << " ";
                std::time(&timer2);
                total_exec_time = std::difftime(timer2, timer1);
                std::cout << "Total computing time: " << total_exec_time << " s; ";
                std::cout << std::endl;
            }
        }
        std::cout << "Total computing time: " << total_exec_time << " s; ";
    } catch (cl::Error & err) {
        std::cout << err.what() << " " << err.err() << std::endl;
    }
}

CLVisc::~CLVisc() {
}

// in case the input array is a vector of double instead of float
template void CLVisc::read_ini(const std::vector<float> & ed);
template void CLVisc::read_ini(const std::vector<double> & ed);

template void CLVisc::read_ini(const std::vector<float> & ed, 
                      const std::vector<float> & vx, 
                      const std::vector<float> & vy, 
                      const std::vector<float> & vz);

template void CLVisc::read_ini(const std::vector<double> & ed, 
                      const std::vector<double> & vx, 
                      const std::vector<double> & vy, 
                      const std::vector<double> & vz);


template void CLVisc::read_ini(const std::vector<float> & ed, 
              const std::vector<float> & vx, 
              const std::vector<float> & vy, 
              const std::vector<float> & vz,
              const std::vector<float> & pi00,
              const std::vector<float> & pi01,
              const std::vector<float> & pi02,
              const std::vector<float> & pi03,
              const std::vector<float> & pi11,
              const std::vector<float> & pi12,
              const std::vector<float> & pi13,
              const std::vector<float> & pi22,
              const std::vector<float> & pi23,
              const std::vector<float> & pi33);

template void CLVisc::read_ini(const std::vector<double> & ed, 
              const std::vector<double> & vx, 
              const std::vector<double> & vy, 
              const std::vector<double> & vz,
              const std::vector<double> & pi00,
              const std::vector<double> & pi01,
              const std::vector<double> & pi02,
              const std::vector<double> & pi03,
              const std::vector<double> & pi11,
              const std::vector<double> & pi12,
              const std::vector<double> & pi13,
              const std::vector<double> & pi22,
              const std::vector<double> & pi23,
              const std::vector<double> & pi33);


} // end namespace clvisc
