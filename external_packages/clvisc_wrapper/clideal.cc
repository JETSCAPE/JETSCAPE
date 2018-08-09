#include <cmath>
#include <cstring>
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

    opts_.KernelIncludePath("../../PyVisc/pyvisc/kernel/");
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
    read_eos_table_("../eos_table/s95_pce165.dat");
    try {
        // build kernels for hydrodynamic evolution
        auto prg = backend_.BuildProgram("../../PyVisc/pyvisc/kernel/kernel_ideal.cl", opts_);
        kernel_kt_src_christoffel_ = cl::Kernel(prg, "kt_src_christoffel");
        kernel_kt_src_alongx_ = cl::Kernel(prg, "kt_src_alongx");
        kernel_kt_src_alongy_ = cl::Kernel(prg, "kt_src_alongy");
        kernel_kt_src_alongz_ = cl::Kernel(prg, "kt_src_alongz");
        kernel_update_ev_ = cl::Kernel(prg, "update_ev");
        // build kernels to look for maximum energy density
        auto prg2 = backend_.BuildProgram("../../PyVisc/pyvisc/kernel/kernel_reduction.cl",
                                          opts_);
        kernel_reduction_ = cl::Kernel(prg2, "reduction_stage1");
    } catch (cl::Error & err ){
        std::cerr<<"Error:"<<err.what()<<"("<<err.err()<<")\n";
    }

    //load_initial_condition(data);
    //create gpu buffer

}


void CLIdeal::read_eos_table_(std::string fname) {
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
            host_eos_table.push_back((cl_float4){speed_of_sound_square,
                    pressure, temperature, entropy_density});
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
}

// read initial energy density from external vector
void CLIdeal::read_ini_ed(const std::vector<cl_real> & ed) {
    for (size_t idx = 0; idx < ed.size(); idx++) {
        h_ev_.push_back((cl_real4){ed.at(idx), 0.0f, 0.0f, 0.0f});
    }
    initialize_gpu_buffer_();
}

// read initial ed, vx, vy, vz vector
void CLIdeal::read_ini(const std::vector<cl_real> & ed, 
                       const std::vector<cl_real> & vx, 
                       const std::vector<cl_real> & vy, 
                       const std::vector<cl_real> & vz)
{
    for (size_t idx = 0; idx < ed.size(); idx++) {
        h_ev_.push_back((cl_real4){ed.at(idx), vx.at(idx), vy.at(idx), vz.at(idx)});
    }
    initialize_gpu_buffer_();
}


// step update for Runge-Kutta method, step = {1, 2}
void CLIdeal::step_update(int step) {
    auto size = cfg_.nx * cfg_.ny * cfg_.nz;
    kernel_kt_src_christoffel_.setArg(0, d_src_);
    kernel_kt_src_christoffel_.setArg(1, d_ev_[step]);
    kernel_kt_src_christoffel_.setArg(2, eos_table_);
    kernel_kt_src_christoffel_.setArg(3, tau_);
    kernel_kt_src_christoffel_.setArg(4, step);
    cl::Event event;
    backend_.Queue().enqueueNDRangeKernel(
            kernel_kt_src_christoffel_,        // kernel name
            cl::NullRange,                     // offset 
            cl::NDRange(size),                 // global size
            cl::NullRange,                     // local size (automatically set by system)
            NULL,                              // event waitting list
            &event);                           // event for profiling
    event.wait();
    auto exec_time = backend_.ExcutionTime(event);
    std::cout << "Excution time for christoffel is: " << exec_time << std::endl;

    // update along x direction
    kernel_kt_src_alongx_.setArg(0, d_src_);
    kernel_kt_src_alongx_.setArg(1, d_ev_[step]);
    kernel_kt_src_alongx_.setArg(2, eos_table_);
    kernel_kt_src_alongx_.setArg(3, tau_);
    cl::Event event_alongx;
    backend_.Queue().enqueueNDRangeKernel(
            kernel_kt_src_alongx_,        // kernel name
            cl::NullRange,                     // offset 
            cl::NDRange(cfg_.block_size, cfg_.ny, cfg_.nz),  // global size
            cl::NDRange(cfg_.block_size, 1, 1), // local size (automatically set by system)
            NULL,                              // event waitting list
            &event_alongx);                           // event for profiling
    event_alongx.wait();
    exec_time = backend_.ExcutionTime(event_alongx);
    std::cout << "Excution time for x update is: " << exec_time << std::endl;

    // update along y direction
    kernel_kt_src_alongy_.setArg(0, d_src_);
    kernel_kt_src_alongy_.setArg(1, d_ev_[step]);
    kernel_kt_src_alongy_.setArg(2, eos_table_);
    kernel_kt_src_alongy_.setArg(3, tau_);
    cl::Event event_alongy;
    backend_.Queue().enqueueNDRangeKernel(
            kernel_kt_src_alongy_,        // kernel name
            cl::NullRange,                     // offset 
            cl::NDRange(cfg_.nx, cfg_.block_size, cfg_.nz),  // global size
            cl::NDRange(1, cfg_.block_size, 1), // local size (automatically set by system)
            NULL,                              // event waitting list
            &event_alongy);                           // event for profiling
    event_alongy.wait();
    exec_time = backend_.ExcutionTime(event_alongy);
    std::cout << "Excution time for y update is: " << exec_time << std::endl;

    // update along spacetime-rapidity direction
    kernel_kt_src_alongz_.setArg(0, d_src_);
    kernel_kt_src_alongz_.setArg(1, d_ev_[step]);
    kernel_kt_src_alongz_.setArg(2, eos_table_);
    kernel_kt_src_alongz_.setArg(3, tau_);
    cl::Event event_alongz;
    backend_.Queue().enqueueNDRangeKernel(
            kernel_kt_src_alongz_,        // kernel name
            cl::NullRange,                     // offset 
            cl::NDRange(cfg_.nx, cfg_.ny, cfg_.block_size),  // global size
            cl::NDRange(1, 1, cfg_.block_size), // local size (automatically set by system)
            NULL,                              // event waitting list
            &event_alongz);                           // event for profiling
    event_alongz.wait();
    exec_time = backend_.ExcutionTime(event_alongz);
    std::cout << "Excution time for z update is: " << exec_time << std::endl;

    // update energy density and fluid velocity for all cells
    kernel_update_ev_.setArg(0, d_ev_[3-step]);
    kernel_update_ev_.setArg(1, d_ev_[1]);
    kernel_update_ev_.setArg(2, d_src_);
    kernel_update_ev_.setArg(3, eos_table_);
    kernel_update_ev_.setArg(4, tau_);
    kernel_update_ev_.setArg(5, step);
    cl::Event event_update_ev;
    backend_.Queue().enqueueNDRangeKernel(
            kernel_update_ev_,        // kernel name
            cl::NullRange,            // offset 
            cl::NDRange(size),        // global size
            cl::NullRange,            // local size (automatically set by system)
            NULL,                     // event waitting list
            &event_update_ev);           // event for profiling
    event_update_ev.wait();
    exec_time = backend_.ExcutionTime(event_update_ev);
    std::cout << "Excution time for update e v: " << exec_time << std::endl;
};


void CLIdeal::evolve() {
    int max_loops = 2000;
    for (int step=0; step < max_loops; step++) {
        std::cout << "tau = " << tau_ << " fm" << std::endl;
        step_update(1);
        step_update(2);
        tau_ += cfg_.dt;
        std::cout << std::endl;
    }
}

CLIdeal::~CLIdeal() {
}

} // end namespace clvisc
