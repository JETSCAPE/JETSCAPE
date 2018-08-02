#include "include/opencl_backend.h"
#include <iostream>

int main(int argc, char ** argv) {
    auto backend = OpenclBackend("cpu", 0);
    return 0;
}
