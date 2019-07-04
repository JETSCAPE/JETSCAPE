#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <vector>
#include "lorentz.h"

double delta_qhat(int pid, double E, double M, double T);
double qhat_small_angle_LOpQCD(int pid, double E, double M, double T);
double qhat_L_small_angle_LOpQCD(int pid, double E, double M, double T);
double qhat(int pid, double E, double M, double T);
double qhat_L(int pid, double E, double M, double T);
double dqhat_L_dp2(int pid, double E, double M, double T);
void Ito_update(int pid, double dt, double M, double T, std::vector<double> v, const fourvec & pIn, fourvec & pOut);

#endif
