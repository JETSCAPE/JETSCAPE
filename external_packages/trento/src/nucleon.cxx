// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// TRENTO3D: Three-dimensional extension of TRENTO by Weiyao Ke
// MIT License

#include "nucleon.h"

#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/program_options/variables_map.hpp>

#include "fwd_decl.h"

namespace trento {

namespace {

// Create ctor parameters for unit mean std::gamma_distribution.
//   mean = alpha*beta == 1  ->  beta = 1/alpha
// Used below in NucleonProfile ctor initializer list.

template <typename RealType> using param_type =
  typename std::gamma_distribution<RealType>::param_type;

template <typename RealType>
param_type<RealType> gamma_param_unit_mean(RealType alpha = 1.) {
  return param_type<RealType>{alpha, 1./alpha};
}

// These constants define distances in terms of the width of the nucleon profile
// Gaussian thickness function.

// Truncation radius of the thickness function.
constexpr double trunc_radius_widths = 5.;

// Maximum impact parameter for participation.
constexpr double max_impact_widths = 6.;

// Trivial helper function.
template <typename T>
constexpr T sqr(T value) {
  return value * value;
}

// Inelastic nucleon-nucleon cross section as function of beam energy sqrt(s)
// Fit coefficients explained in the docs.
double cross_sec_from_energy(double sqrts) {
  auto a = 3.1253;
  auto b = 0.1280;
  auto c = 2.0412;
  auto d = 1.8231;
  return a + b * pow(std::log(sqrts) - c, d);
}

// Determine the cross section parameter for sampling participants.
// See section "Fitting the cross section" in the online docs.
double compute_cross_sec_param(const VarMap& var_map) {
  // Read parameters from the configuration.

  // Use manual inelastic nucleon-nucleon cross section if specified.
  // Otherwise default to extrapolated cross section.
  auto sigma_nn = var_map["cross-section"].as<double>();
  if (sigma_nn < 0) {
    sigma_nn = cross_sec_from_energy(var_map["beam-energy"].as<double>());
  }
  auto width = var_map["nucleon-width"].as<double>();

  // Initialize arguments for boost root finding function.

  // Bracket min and max.
  auto a = -10.;
  auto b = 20.;

  // Tolerance function.
  // Require 3/4 of double precision.
  math::tools::eps_tolerance<double> tol{
    (std::numeric_limits<double>::digits * 3) / 4};

  // Maximum iterations.
  // This is overkill -- in testing only 10-20 iterations were required
  // (but no harm in overestimating).
  boost::uintmax_t max_iter = 1000;

  // The right-hand side of the equation.
  auto rhs = sigma_nn / (4 * math::double_constants::pi * sqr(width));

  // This quantity appears a couple times in the equation.
  auto c = sqr(max_impact_widths) / 4;

  try {
    auto result = math::tools::toms748_solve(
      [&rhs, &c](double x) {
        using std::exp;
        using math::expint;
        return c - expint(-exp(x)) + expint(-exp(x-c)) - rhs;
      },
      a, b, tol, max_iter);

    return .5*(result.first + result.second);
  }
  catch (const std::domain_error&) {
    // Root finding fails for very small nucleon widths, w^2/sigma_nn < ~0.01.
    throw std::domain_error{
      "unable to fit cross section -- nucleon width too small?"};
  }
}

}  // unnamed namespace

NucleonProfile::NucleonProfile(const VarMap& var_map)
    : width_sqr_(sqr(var_map["nucleon-width"].as<double>())),
      trunc_radius_sqr_(sqr(trunc_radius_widths)*width_sqr_),
      max_impact_sqr_(sqr(max_impact_widths)*width_sqr_),
      neg_one_div_two_width_sqr_(-.5/width_sqr_),
	  neg_one_div_four_width_sqr_(-.25/width_sqr_),
	  one_div_four_pi_(0.5*math::double_constants::one_div_two_pi),
      cross_sec_param_(compute_cross_sec_param(var_map)),
      fast_exp_(-.5*sqr(trunc_radius_widths), 0., 1000),
      fluct_dist_(gamma_param_unit_mean(var_map["fluctuation"].as<double>())),
      prefactor_(math::double_constants::one_div_two_pi/width_sqr_),
      with_ncoll_(var_map["ncoll"].as<bool>())
{}

}  // namespace trento
