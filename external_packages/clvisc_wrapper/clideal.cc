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
    read_eos_table("../eos_table/s95_pce165.dat");
    auto prg = backend_.BuildProgram("../../PyVisc/pyvisc/kernel/kernel_ideal.cl", opts_);
    try {
        kernel_kt_src_christofeel_ = cl::Kernel(prg, "kt_src_christoffel");
        kernel_kt_src_alongx_ = cl::Kernel(prg, "kt_src_alongx");
        kernel_kt_src_alongy_ = cl::Kernel(prg, "kt_src_alongy");
        kernel_kt_src_alongz_ = cl::Kernel(prg, "kt_src_alongz");
        kernel_update_ev_ = cl::Kernel(prg, "update_ev");
    } catch (cl::Error & err ){
        std::cerr<<"Error:"<<err.what()<<"("<<err.err()<<")\n";
    }
}


void CLIdeal::read_eos_table(std::string fname) {
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


CLIdeal::~CLIdeal()
{
}

void CLIdeal::clean()
{
};

