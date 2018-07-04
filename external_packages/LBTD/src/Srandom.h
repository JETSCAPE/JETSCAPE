#ifndef SRANDOM_H
#define SRANDOM_H
#include <random>

namespace Srandom{
extern std::mt19937 gen;
extern std::uniform_real_distribution<double> sqrtZ;
extern std::uniform_real_distribution<double> rejection;
extern std::uniform_real_distribution<double> init_dis;
extern std::uniform_real_distribution<double> dist_phi;
extern std::uniform_real_distribution<double> dist_costheta;
}
#endif
