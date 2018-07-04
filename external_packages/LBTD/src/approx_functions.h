#ifndef APPROX_FUNCTIONS_H
#define APPROX_FUNCTIONS_H
#include <vector>
#include "lorentz.h"

// Xsection
scalar approx_X22(std::vector<double> params);
scalar approx_dX22_max(std::vector<double> params);
scalar approx_X23(std::vector<double> params);
scalar approx_dX23_max(std::vector<double> params);
scalar approx_X32(std::vector<double> params);
scalar approx_dX32_max(std::vector<double> params);
// Rate
scalar approx_R22(std::vector<double> params);
scalar approx_dR22_max(std::vector<double> params);
scalar approx_R23(std::vector<double> params);
scalar approx_dR23_max(std::vector<double> params);
scalar approx_R32(std::vector<double> params);
scalar approx_dR32_max(std::vector<double> params);
#endif
