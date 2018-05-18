#include "random.h"

namespace Srandom{
const double AMC = 4.0;
std::mt19937 gen(std::random_device{}());
std::uniform_real_distribution<double> sqrtZ(std::sqrt(1./AMC), std::sqrt(AMC));
std::uniform_real_distribution<double> rejection(0.0, 1.0);
std::uniform_real_distribution<double> init_dis(0.0, 1.0);
std::uniform_real_distribution<double> dist_phi(0.0, 2.0*M_PI);
std::uniform_real_distribution<double> dist_costheta(-1.0, 1.0);
std::normal_distribution<double> white_noise(0.0, 1.0);
}
