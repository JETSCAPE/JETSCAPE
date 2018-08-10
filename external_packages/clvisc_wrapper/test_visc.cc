#include "include/opencl_backend.h"
#include "include/clideal.h"
#include "include/clvisc.h"
#include <iostream>

clvisc::Config set_hydro_config() {
    clvisc::Config cfg;
    cfg.block_size = 64;
    cfg.nx = 200;
    cfg.ny = 200;
    cfg.nz = 1;
    cfg.tau0 = 0.6;
    cfg.dt = 0.02;
    cfg.dx = 0.16;
    cfg.dy = 0.16;
    cfg.dz = 0.16;
    cfg.etaos_xmin = 0.15;
    cfg.etaos_ymin = 0.15;
    cfg.etaos_left_slop = 0.0;
    cfg.etaos_right_slop = 0.0;
    cfg.result_directory = "./";
    return cfg;
}


std::vector<clvisc::cl_real> trivial_initial_condition(const clvisc::Config & cfg) {
    std::vector<clvisc::cl_real> ini;
    long ngrids = cfg.nx * cfg.ny * cfg.nz;
    ini.reserve(ngrids);
    for (long idx = 0; idx < ngrids; idx++) {
        clvisc::cl_real const_ed = 30.0f;
        ini.push_back(const_ed);
    }
    return ini;
}

int main(int argc, char ** argv) {
    auto backend0 = clvisc::OpenclBackend("cpu", 0);
    backend0.DeviceInfo();

    auto backend1 = clvisc::OpenclBackend("gpu", 1);
    backend1.DeviceInfo();

    auto cfg = set_hydro_config();
    std::string device_type = "gpu";
    int device_id = 1;
    clvisc::CLVisc visc(cfg, device_type, device_id);
    auto ini = trivial_initial_condition(cfg);
    visc.read_ini(ini);
    visc.evolve();
    return 0;
}
