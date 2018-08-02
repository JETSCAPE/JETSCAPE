#include "include/opencl_backend.h"
#include "include/clideal.h"
#include <iostream>

int main(int argc, char ** argv) {
    auto backend0 = OpenclBackend("cpu", 0);
    backend0.DeviceInfo();

    auto backend1 = OpenclBackend("gpu", 1);
    backend1.DeviceInfo();

    CLIdeal ideal("./", "gpu", 0);
    return 0;
}
