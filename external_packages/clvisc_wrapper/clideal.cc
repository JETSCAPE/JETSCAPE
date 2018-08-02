#include <cmath>
#include <cstring>
#include "include/clideal.h"
#define REDUCTION_BLOCKS 64


CLIdeal::CLIdeal(std::string config_file_path, std::string device_type,
        int device_id):backend_(OpenclBackend(device_type, device_id)), 
                       data_path_(config_file_path){
    int NX = 200;
    int NY = 200;
    int NZ = 200;
    int SIZE = NX * NY * NZ;
    cl_real tau0 = 0.6;
    cl_real dt = 0.02;
    cl_real dx = 0.2;
    cl_real dy = 0.2;
    cl_real dz = 0.2;
    cl_real etaos_xmin = 0.15;
    cl_real etaos_ymin = 0.15;
    cl_real etaos_left_slop = 0.0;
    cl_real etaos_right_slop = 0.0;
    cl_real lambda1 = -10.0;
    opts_.KernelIncludePath("../../PyVisc/pyvisc/kernel/");
    opts_.Define("EOS_TABLE");
    opts_.SetFloatConst("TAU0", tau0);
    opts_.SetFloatConst("DT", dt);
    opts_.SetFloatConst("DX", dx);
    opts_.SetFloatConst("DY", dy);
    opts_.SetFloatConst("DZ", dz);
    opts_.SetFloatConst("ETAOS_XMIN", etaos_xmin);
    opts_.SetFloatConst("ETAOS_YMIN", etaos_ymin);
    opts_.SetFloatConst("ETAOS_LEFT_SLOP", etaos_left_slop);
    opts_.SetFloatConst("ETAOS_RIGHT_SLOP", etaos_right_slop);
    opts_.SetFloatConst("LAM1", lambda1);
    opts_.SetIntConst("BSZ", BSZ);
    opts_.SetIntConst("NX", NX);
    opts_.SetIntConst("NY", NY);
    opts_.SetIntConst("NZ", NZ);
    opts_.SetIntConst("SIZE", SIZE);
    auto prg = backend_.BuildProgram("../../PyVisc/pyvisc/kernel/kernel_ideal.cl", opts_);
    kernel_kt_src_christofeel_ = cl::Kernel(prg, "kernel_kt_src_christofeel");
    kernel_kt_src_alongx_ = cl::Kernel(prg, "kernel_kt_src_alongx");
    kernel_kt_src_alongy_ = cl::Kernel(prg, "kernel_kt_src_alongy");
    kernel_kt_src_alongz_ = cl::Kernel(prg, "kernel_kt_src_alongz");
    kernel_update_ev_ = cl::Kernel(prg, "kernel_update_ev_");
}


CLIdeal::~CLIdeal()
{
}

void CLIdeal::clean()
{
};

