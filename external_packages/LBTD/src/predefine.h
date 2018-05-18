#ifndef PREDEFINE_H
#define PREDEFINE_H

#include <cmath>
#include <H5Cpp.h>

//=============constants=======================================================
const double c4d9 = 4./9.;
const double c1d9 = 1./9.;
const double c16pi = 16.*M_PI;
const double c48pi = 48.*M_PI;
const double c16pi2 = 16.*M_PI*M_PI;
const double c64d9pi2 = 64./9.*M_PI*M_PI;
const double c256pi4 = 256.*std::pow(M_PI, 4);
const int Nc = 3, nf = 3;
const double pf_g = 4.*M_PI/3.*(Nc + nf/2.); // prefractor for gluon self energy^2
const double pf_q = M_PI/2.*(Nc*Nc - 1)/2./Nc; // prefractor for quark self energy^2
const double alpha0 = 4.*M_PI/(11. - 2./3.*nf); // alpha_s(Q2 = e*Lambda2)
const double Lambda = 0.2; // [GeV] Lambda QCD
const double Lambda2 = Lambda*Lambda; // [GeV^2] Lambda QCD squared
const double mu2_left = Lambda2*std::exp(1.1),
			 mu2_right = Lambda2*1e-2;
const double fmc_to_GeV_m1 = 5.026;

// [GeV^2] ranges within which alphas > 1 and will be cut

template <typename T> inline const H5::PredType& type();
template <> inline const H5::PredType& type<size_t>() { return H5::PredType::NATIVE_HSIZE; }
template <> inline const H5::PredType& type<double>() { return H5::PredType::NATIVE_DOUBLE; }

template <typename T>
void hdf5_add_scalar_attr(
  const H5::Group& gp, const std::string& name, const T& value) {
  const auto& datatype = type<T>();
  auto attr = gp.createAttribute(name.c_str(), datatype, H5::DataSpace{});
  attr.write(datatype, &value);
}

template <typename T>
void hdf5_read_scalar_attr(
  const H5::Group& gp, const std::string& name, T& value) {
  const auto& datatype = type<T>();
  auto attr = gp.openAttribute(name.c_str());
  attr.read(datatype, &value);
}
#endif
