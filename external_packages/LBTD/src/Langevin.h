#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <vector>
#include "lorentz.h"

extern double A, B;
double kperp(double E, double M, double T);
double kpara(double E, double M, double T);
void initialize_transport_coeff(double A, double B);
void postpoint_update( double dt, double M, double T, std::vector<double> v, const fourvec & pIn, fourvec & pOut);
void Ito_update( double dt, double M, double T, std::vector<double> v, const fourvec & pIn, fourvec & pOut);
#endif
