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
#include "include/clideal.h"

namespace clvisc {
// number of blocks for reduction on gpus
const int REDUCTION_BLOCKS = 64;

CLIdeal::CLIdeal(const Config & cfg, std::string device_type,
        int device_id):cfg_(cfg),
        backend_(OpenclBackend(device_type, device_id))
{
    // set current time to initial time
    tau_ = cfg.tau0;

    CompileOption opts_;
    opts_.KernelIncludePath("clvisc_kernel/");
    //opts_.Define("EOSI");
    opts_.Define("EOS_TABLE");
    opts_.SetFloatConst("TAU0", cfg_.tau0);
    opts_.SetFloatConst("DT", cfg_.dt);
    opts_.SetFloatConst("DX", cfg_.dx);
    opts_.SetFloatConst("DY", cfg_.dy);
    opts_.SetFloatConst("DZ", cfg_.dz);
    opts_.SetFloatConst("ETAOS_XMIN", cfg_.etaos_xmin);
    opts_.SetFloatConst("ETAOS_YMIN", cfg_.etaos_ymin);
    opts_.SetFloatConst("ETAOS_LEFT_SLOP", cfg_.etaos_left_slop);
    opts_.SetFloatConst("ETAOS_RIGHT_SLOP", cfg_.etaos_right_slop);
    cl_real lambda1 = -10.0; // one parameters used for analytical solution checking
    opts_.SetFloatConst("LAM1", lambda1);
    opts_.SetIntConst("NX", cfg_.nx);
    opts_.SetIntConst("NY", cfg_.ny);
    opts_.SetIntConst("NZ", cfg_.nz);
    opts_.SetIntConst("BSZ", cfg.block_size);
    int SIZE = cfg_.nx * cfg_.ny * cfg_.nz;
    opts_.SetIntConst("SIZE", SIZE);
#ifdef USE_SINGLE_PRECISION
    opts_.Define("USE_SINGLE_PRECISION");
#endif
    read_eos_table_("data_table/s95_pce165.dat", opts_);
    compile_option_ = opts_.str();
    try {
        // build kernels for hydrodynamic evolution
        auto prg = backend_.BuildProgram("clvisc_kernel/kernel_ideal.cl", compile_option_);
        kernel_kt_src_christoffel_ = cl::Kernel(prg, "kt_src_christoffel");
        kernel_kt_src_alongx_ = cl::Kernel(prg, "kt_src_alongx");
        kernel_kt_src_alongy_ = cl::Kernel(prg, "kt_src_alongy");
        kernel_kt_src_alongz_ = cl::Kernel(prg, "kt_src_alongz");
        kernel_update_ev_ = cl::Kernel(prg, "update_ev");
        // build kernels to look for maximum energy density
        auto prg2 = backend_.BuildProgram("clvisc_kernel/kernel_reduction.cl",
                                          compile_option_);
        kernel_reduction_ = cl::Kernel(prg2, "reduction_stage1");
    } catch (cl::Error & err ){
        std::cerr<<"Error:"<<err.what()<<"("<<err.err()<<")\n";
    }
}


void CLIdeal::read_eos_table_(std::string fname, CompileOption & opts_) {
    // eos_table stored in fname has the following format:
    // cs2, ed [GeV/fm^3], pr[GeV/fm^3], T[GeV]
    std::vector<cl_float4> host_eos_table;
    std::ifstream fin(fname);
    cl_float speed_of_sound_square;
    cl_float pressure;
    cl_float temperature;
    cl_float entropy_density;
    char comments[256];
    if (fin.is_open()) {
        fin.getline(comments, 256);
        while (fin.good()) {
            fin >> speed_of_sound_square >> pressure >> temperature >> entropy_density;
            if (fin.eof()) break;  // eof() repeat the last line
            host_eos_table.push_back((cl_float4){{speed_of_sound_square,
                    pressure, temperature, entropy_density}});
        }
    } else {
        throw std::runtime_error("Failed to open equation of state table" + fname);
    }

    bool read_only = true;
    size_t width = 1555;
    size_t height = 100;

    eos_table_ = backend_.CreateImage2DByCopyVector(host_eos_table,
                          width, height, read_only);
    opts_.SetFloatConst("EOS_ED_START", 0.0);
    opts_.SetFloatConst("EOS_ED_STEP", 0.002);
    opts_.SetIntConst("EOS_NUM_ED", 155500);
    opts_.SetIntConst("EOS_NUM_OF_ROWS", height);
    opts_.SetIntConst("EOS_NUM_OF_COLS", width);
}

void CLIdeal::initialize_gpu_buffer_() {
    bool read_only_option = false;
    for (int i = 0; i < 3; i++) {
        d_ev_[i] = backend_.CreateBufferByCopyVector(h_ev_, read_only_option);
    }
    d_src_ = backend_.CreateBufferByCopyVector(h_ev_, read_only_option);

    d_submax_ = backend_.CreateBuffer(REDUCTION_BLOCKS*sizeof(cl_real));
}

// read initial energy density from external vector
template <typename ValueType>
void CLIdeal::read_ini(const std::vector<ValueType> & ed) {
    h_ev_.clear();
    for (size_t idx = 0; idx < ed.size(); idx++) {
        h_ev_.push_back((cl_real4){{static_cast<cl_real>(ed.at(idx)), 0.0f, 0.0f, 0.0f}});
    }
    initialize_gpu_buffer_();
}

// read initial ed, vx, vy, vz vector
template <typename ValueType>
void CLIdeal::read_ini(const std::vector<ValueType> & ed, 
                       const std::vector<ValueType> & vx, 
                       const std::vector<ValueType> & vy, 
                       const std::vector<ValueType> & vz)
{
    h_ev_.clear();
    for (size_t idx = 0; idx < ed.size(); idx++) {
        h_ev_.push_back((cl_real4){{static_cast<cl_real>(ed.at(idx)), \
                                   static_cast<cl_real>(vx.at(idx)), \
                                   static_cast<cl_real>(vy.at(idx)), \
                                   static_cast<cl_real>(vz.at(idx))}});
    }
    initialize_gpu_buffer_();
}


// step update for Runge-Kutta method, step = {1, 2}
void CLIdeal::half_step_(int step) {
    auto size = cfg_.nx * cfg_.ny * cfg_.nz;
    kernel_kt_src_christoffel_.setArg(0, d_src_);
    kernel_kt_src_christoffel_.setArg(1, d_ev_[step]);
    kernel_kt_src_christoffel_.setArg(2, eos_table_);
    // notice: it is important to cast tau_ from double to float explicitly
    // otherwise the program will produce wrong results
    kernel_kt_src_christoffel_.setArg(3, static_cast<cl_real>(tau_));
    kernel_kt_src_christoffel_.setArg(4, step);
    backend_.enqueue_run(kernel_kt_src_christoffel_,
            cl::NDRange(size),        // global size
            cl::NullRange);           // local size
 
    // update along x direction
    kernel_kt_src_alongx_.setArg(0, d_src_);
    kernel_kt_src_alongx_.setArg(1, d_ev_[step]);
    kernel_kt_src_alongx_.setArg(2, eos_table_);
    kernel_kt_src_alongx_.setArg(3, static_cast<cl_real>(tau_));
    backend_.enqueue_run(kernel_kt_src_alongx_,
            cl::NDRange(cfg_.block_size, cfg_.ny, cfg_.nz),  // global size
            cl::NDRange(cfg_.block_size, 1, 1));             // local size

    // update along y direction
    kernel_kt_src_alongy_.setArg(0, d_src_);
    kernel_kt_src_alongy_.setArg(1, d_ev_[step]);
    kernel_kt_src_alongy_.setArg(2, eos_table_);
    kernel_kt_src_alongy_.setArg(3, static_cast<cl_real>(tau_));
    backend_.enqueue_run(kernel_kt_src_alongy_,
            cl::NDRange(cfg_.nx, cfg_.block_size, cfg_.nz),  // global size
            cl::NDRange(1, cfg_.block_size, 1));             // local size

    if (cfg_.nz != 1) {
        // update along spacetime-rapidity direction
        kernel_kt_src_alongz_.setArg(0, d_src_);
        kernel_kt_src_alongz_.setArg(1, d_ev_[step]);
        kernel_kt_src_alongz_.setArg(2, eos_table_);
        kernel_kt_src_alongz_.setArg(3, static_cast<cl_real>(tau_));
        backend_.enqueue_run(kernel_kt_src_alongz_,
                cl::NDRange(cfg_.nx, cfg_.ny, cfg_.block_size),  // global size
                cl::NDRange(1, 1, cfg_.block_size));             // local size
    }

    //// update energy density and fluid velocity for all cells
    kernel_update_ev_.setArg(0, d_ev_[3-step]);
    kernel_update_ev_.setArg(1, d_ev_[1]);
    kernel_update_ev_.setArg(2, d_src_);
    kernel_update_ev_.setArg(3, eos_table_);
    kernel_update_ev_.setArg(4, static_cast<cl_real>(tau_));
    kernel_update_ev_.setArg(5, step);
    backend_.enqueue_run(kernel_update_ev_,
            cl::NDRange(size),        // global size
            cl::NullRange);           // local size
            
    //backend_.enqueue_copy(d_src_, h_ev_);
    //std::cout << h_ev_.at(0).s[0] << " ";
    //std::cout << std::endl;
}

// return the maximum energy density of the QGP
// will be used to stop hydro
float CLIdeal::max_energy_density() {
    int size = cfg_.nx * cfg_.ny * cfg_.nz;
    kernel_reduction_.setArg(0, d_ev_[1]);
    kernel_reduction_.setArg(1, d_submax_);
    kernel_reduction_.setArg(2, size);
    auto global_size = cl::NDRange(256*REDUCTION_BLOCKS);
    auto local_size = cl::NDRange(256);
    backend_.enqueue_run(kernel_reduction_, global_size, local_size);

    std::vector<cl_real> h_submax(REDUCTION_BLOCKS);
    backend_.enqueue_copy(d_submax_, h_submax);
    return *std::max_element(h_submax.begin(), h_submax.end());
}

// run hydrodynamic evolution for one time step
// return the excution time on device
void CLIdeal::one_step() {
    half_step_(1);
    half_step_(2);
}

// predict the first half step; usful to get initial fluid velocity
// for viscous hydrodynamics 
void CLIdeal::predict_first_step() {
    half_step_(1);
}

void CLIdeal::evolve() {
    int max_loops = 2000;
    float total_exec_time = 0.0;
    try {
        double total_exec_time = 0.0;
        std::time_t timer1, timer2;
        std::time(&timer1);
        for (int step=0; step < max_loops; step++) {
            float max_ed = max_energy_density();
            if ( max_ed < 0.5 ) break;
            one_step();
            std::time(&timer2);
            float time_diff = std::difftime(timer2, timer1);
            std::cout << "tau = " << tau_ << " fm; ";
            std::cout << "max_ed = " << max_ed << " ";
            std::cout << "Total computing time: " << time_diff << " s; ";
            std::cout << std::endl;
            tau_ += cfg_.dt;
        }
        std::cout << "Total computing time: " << total_exec_time << " s; ";
    } catch (cl::Error & err) {
        std::cout << err.what() << " " << err.err() << std::endl;
    }
}

OpenclBackend & CLIdeal::get_backend() {
    return backend_;
}

std::string & CLIdeal::get_compile_option() {
    return compile_option_;
}


CLIdeal::~CLIdeal() {
}



template void CLIdeal::read_ini(const std::vector<float> & ed);
template void CLIdeal::read_ini(const std::vector<double> & ed);

template void CLIdeal::read_ini(const std::vector<float> & ed, 
                       const std::vector<float> & vx, 
                       const std::vector<float> & vy, 
                       const std::vector<float> & vz);

template void CLIdeal::read_ini(const std::vector<double> & ed, 
                       const std::vector<double> & vx, 
                       const std::vector<double> & vy, 
                       const std::vector<double> & vz);


} // end namespace clvisc
