// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "nucleus.h"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <boost/math/constants/constants.hpp>
#ifdef TRENTO_HDF5
// include multi_array for use with ManualNucleus
#ifdef NDEBUG
#define BOOST_DISABLE_ASSERTS
#endif
#include <boost/multi_array.hpp>
#endif

#include "hdf5_utils.h"
#include "random.h"

namespace trento {

namespace {

// Correct Woods-Saxon surface thickness parameter (a) for finite Gaussian
// nucleon width (w):
//
//    a_corrected^2 = a^2 - c^2*w*2
//
// where c is a universal constant independent of a and w.
//
// See https://gist.github.com/jbernhard/60b3ab9662a4737658d8.
double correct_a(double a, double w) {
  constexpr auto c = 0.61;  // correction coefficient
  constexpr auto a_min = 0.01;  // min. value (prevent div. by zero, etc.)
  return std::sqrt(std::fmax(a*a - c*c*w*w, a_min*a_min));
}

}  // unnamed namespace

NucleusPtr Nucleus::create(const std::string& species, double nucleon_width) {
  // W-S params ref. in header
  // XXX: remember to add new species to the help output in main() and the readme
  if (species == "p")
    return NucleusPtr{new Proton{}};
  else if (species == "d")
    return NucleusPtr{new Deuteron{}};
  else if (species == "Cu")
    return NucleusPtr{new WoodsSaxonNucleus{
       62, 4.20, correct_a(0.596, nucleon_width)
    }};
  else if (species == "Cu2")
    return NucleusPtr{new DeformedWoodsSaxonNucleus{
       62, 4.20, correct_a(0.596, nucleon_width), 0.162, -0.006
    }};
  else if (species == "Au")
    return NucleusPtr{new WoodsSaxonNucleus{
      197, 6.38, correct_a(0.535, nucleon_width)
    }};
  else if (species == "Au2")
    return NucleusPtr{new DeformedWoodsSaxonNucleus{
      197, 6.38, correct_a(0.535, nucleon_width), -0.131, -0.031
    }};
  else if (species == "Pb")
    return NucleusPtr{new WoodsSaxonNucleus{
      208, 6.62, correct_a(0.546, nucleon_width)
    }};
  else if (species == "U")
    return NucleusPtr{new DeformedWoodsSaxonNucleus{
      238, 6.81, correct_a(0.600, nucleon_width), 0.280, 0.093
    }};
  else if (species == "U2")
    return NucleusPtr{new DeformedWoodsSaxonNucleus{
      238, 6.86, correct_a(0.420, nucleon_width), 0.265, 0.000
    }};
  else if (species == "U3")
    return NucleusPtr{new DeformedWoodsSaxonNucleus{
      238, 6.67, correct_a(0.440, nucleon_width), 0.280, 0.093
    }};
  // Read nuclear configurations from HDF5.
  else if (hdf5::filename_is_hdf5(species)) {
#ifdef TRENTO_HDF5
    return ManualNucleus::create(species);
#else
    throw std::invalid_argument{"HDF5 output was not compiled"};
#endif  // TRENTO_HDF5
  }
  else
    throw std::invalid_argument{"unknown projectile species: " + species};
}

Nucleus::Nucleus(std::size_t A) : nucleons_(A), offset_(0) {}

void Nucleus::sample_nucleons(double offset) {
  offset_ = offset;
  sample_nucleons_impl();
}

void Nucleus::set_nucleon_position(Nucleon& nucleon, double x, double y) {
  nucleon.set_position(x + offset_, y);
}

Proton::Proton() : Nucleus(1) {}

/// Always zero.
double Proton::radius() const {
  return 0.;
}

/// Always place the nucleon at the origin.
void Proton::sample_nucleons_impl() {
  set_nucleon_position(*begin(), 0., 0.);
}

// Without loss of generality, let the internal a_ parameter be the minimum of
// the given (a, b) and the internal b_ be the maximum.
Deuteron::Deuteron(double a, double b)
    : Nucleus(2),
      a_(std::fmin(a, b)),
      b_(std::fmax(a, b))
{}

double Deuteron::radius() const {
  // The quantile function for the exponential distribution exp(-2*a*r) is
  // -log(1-q)/(2a).  Return the 99% quantile.
  return -std::log(.01)/(2*a_);
}

void Deuteron::sample_nucleons_impl() {
  // Sample the inter-nucleon radius using rejection sampling with an envelope
  // function.  The Hulthén wavefunction including the r^2 Jacobian expands to
  // three exponential terms:  exp(-2*a*r) + exp(-2*b*r) - 2*exp(-(a+b)*r).
  // This does not have a closed-form inverse CDF, however we can easily sample
  // exponential numbers from the term that falls off the slowest, i.e.
  // exp(-2*min(a,b)*r).  In the ctor initializer list the "a" parameter is
  // always set to the minimum, so we should sample from exp(-2*a*r).
  double r, prob;
  do {
    // Sample a uniform random number, u = exp(-2*a*r).
    auto u = random::canonical<double>();
    // Invert to find the actual radius.
    r = -std::log(u) / (2*a_);
    // The acceptance probability is now the radial wavefunction over the
    // envelope function, both evaluated at the proposal radius r.
    // Conveniently, the envelope evaluated at r is just the uniform random
    // number u.
    prob = std::pow(std::exp(-a_*r) - std::exp(-b_*r), 2) / u;
  } while (prob < random::canonical<double>());

  // Now sample spherical rotation angles.
  auto cos_theta = random::cos_theta<double>();
  auto phi = random::phi<double>();

  // And compute the transverse coordinates of one nucleon.
  auto r_sin_theta = r * std::sqrt(1. - cos_theta*cos_theta);
  auto x = r_sin_theta * std::cos(phi);
  auto y = r_sin_theta * std::sin(phi);

  // Place the first nucleon at the sampled coordinates (x, y).
  set_nucleon_position(*begin(), x, y);
  // Place the second nucleon opposite to the first, at (-x, -y).
  set_nucleon_position(*std::next(begin()), -x, -y);
}

// Extend the W-S dist out to R + 10a; for typical values of (R, a), the
// probability of sampling a nucleon beyond this radius is O(10^-5).
WoodsSaxonNucleus::WoodsSaxonNucleus(std::size_t A, double R, double a)
    : Nucleus(A),
      R_(R),
      a_(a),
      woods_saxon_dist_(1000, 0., R + 10.*a,
        [R, a](double r) { return r*r/(1.+std::exp((r-R)/a)); })
{}

/// Return something a bit smaller than the true maximum radius.  The
/// Woods-Saxon distribution falls off very rapidly (exponentially), and since
/// this radius determines the impact parameter range, the true maximum radius
/// would cause far too many events with zero participants.
double WoodsSaxonNucleus::radius() const {
  return R_ + 3.*a_;
}

/// Sample uncorrelated Woods-Saxon nucleon positions.
void WoodsSaxonNucleus::sample_nucleons_impl() {
  for (auto&& nucleon : *this) {
    // Sample spherical radius from Woods-Saxon distribution.
    auto r = woods_saxon_dist_(random::engine);

    // Sample isotropic spherical angles.
    auto cos_theta = random::cos_theta<double>();
    auto phi = random::phi<double>();

    // Convert to transverse Cartesian coordinates
    auto r_sin_theta = r * std::sqrt(1. - cos_theta*cos_theta);
    auto x = r_sin_theta * std::cos(phi);
    auto y = r_sin_theta * std::sin(phi);

    set_nucleon_position(nucleon, x, y);
  }
  // XXX: re-center nucleon positions?
}

// Set rmax like the non-deformed case (R + 10a), but for the maximum
// "effective" radius.  The numerical coefficients for beta2 and beta4 are the
// approximate values of Y20 and Y40 at theta = 0.
DeformedWoodsSaxonNucleus::DeformedWoodsSaxonNucleus(
    std::size_t A, double R, double a, double beta2, double beta4)
    : Nucleus(A),
      R_(R),
      a_(a),
      beta2_(beta2),
      beta4_(beta4),
      rmax_(R*(1. + .63*std::fabs(beta2) + .85*std::fabs(beta4)) + 10.*a)
{}

/// Return something a bit smaller than the true maximum radius.  The
/// Woods-Saxon distribution falls off very rapidly (exponentially), and since
/// this radius determines the impact parameter range, the true maximum radius
/// would cause far too many events with zero participants.
double DeformedWoodsSaxonNucleus::radius() const {
  return rmax_ - 7.*a_;
}

double DeformedWoodsSaxonNucleus::deformed_woods_saxon_dist(
    double r, double cos_theta) const {
  auto cos_theta_sq = cos_theta*cos_theta;

  // spherical harmonics
  using math::double_constants::one_div_root_pi;
  auto Y20 = std::sqrt(5)/4. * one_div_root_pi * (3.*cos_theta_sq - 1.);
  auto Y40 = 3./16. * one_div_root_pi *
             (35.*cos_theta_sq*cos_theta_sq - 30.*cos_theta_sq + 3.);

  // "effective" radius
  auto Reff = R_ * (1. + beta2_*Y20 + beta4_*Y40);

  return 1. / (1. + std::exp((r - Reff) / a_));
}

/// Sample uncorrelated deformed Woods-Saxon nucleon positions.
void DeformedWoodsSaxonNucleus::sample_nucleons_impl() {
  // The deformed W-S distribution is defined so the symmetry axis is aligned
  // with the Z axis, so e.g. the long axis of uranium coincides with Z.
  //
  // After sampling positions, they must be randomly rotated.  In general this
  // requires three Euler rotations, but in this case we only need two
  // because there is no use in rotating about the nuclear symmetry axis.
  //
  // The two rotations are:
  //  - a polar "tilt", i.e. rotation about the X axis
  //  - an azimuthal "spin", i.e. rotation about the original Z axis

  // "tilt" angle
  const auto cos_a = random::cos_theta<double>();
  const auto sin_a = std::sqrt(1. - cos_a*cos_a);

  // "spin" angle
  const auto angle_b = random::phi<double>();
  const auto cos_b = std::cos(angle_b);
  const auto sin_b = std::sin(angle_b);

  for (auto&& nucleon : *this) {
    // Sample (r, theta) using a standard rejection method.
    // Remember to include the phase-space factors.
    double r, cos_theta;
    do {
      r = rmax_ * std::cbrt(random::canonical<double>());
      cos_theta = random::cos_theta<double>();
    } while (random::canonical<double>() > deformed_woods_saxon_dist(r, cos_theta));

    // Sample azimuthal angle.
    auto phi = random::phi<double>();

    // Convert to Cartesian coordinates.
    auto r_sin_theta = r * std::sqrt(1. - cos_theta*cos_theta);
    auto x = r_sin_theta * std::cos(phi);
    auto y = r_sin_theta * std::sin(phi);
    auto z = r * cos_theta;

    // Rotate.
    // The rotation formula was derived by composing the "tilt" and "spin"
    // rotations described above.
    auto x_rot = x*cos_b - y*cos_a*sin_b + z*sin_a*sin_b;
    auto y_rot = x*sin_b + y*cos_a*cos_b - z*sin_a*cos_b;

    set_nucleon_position(nucleon, x_rot, y_rot);
  }
}

#ifdef TRENTO_HDF5

namespace {

// Read a slice of an HDF5 dataset into a boost::multi_array.
template <typename T, std::size_t FileDims, std::size_t MemDims>
boost::multi_array<T, MemDims>
read_dataset(
    const H5::DataSet& dataset,
    const std::array<hsize_t, FileDims>& count,
    const std::array<hsize_t, FileDims>& start,
    const std::array<hsize_t, MemDims>& shape) {
  boost::multi_array<T, MemDims> array{shape};
  auto filespace = dataset.getSpace();
  filespace.selectHyperslab(H5S_SELECT_SET, count.data(), start.data());
  auto memspace = hdf5::make_dataspace(shape);

  dataset.read(array.data(), hdf5::type<T>(), memspace, filespace);

  return array;
}

}  // unnamed namespace

std::unique_ptr<ManualNucleus> ManualNucleus::create(const std::string& path) {
  auto file = hdf5::try_open_file(path);

  // Check that there is a single dataset in the file.
  // Might relax this constraint in the future.
  if (file.getNumObjs() != 1)
    throw std::invalid_argument{
      "file '" + path + "' must contain exactly one object"
    };

  auto name = file.getObjnameByIdx(0);
#if H5_VERSION_GE(1, 8, 13)
  if (file.childObjType(name) != H5O_TYPE_DATASET)  // added v1.8.13
#else
  if (file.getObjTypeByIdx(0) != H5G_DATASET)  // deprecated fall back
#endif
    throw std::invalid_argument{
      "object '" + name + "' in file '" + path + "' is not a dataset"
    };

  // Make dataset object in a unique_ptr for eventual passing to ctor.
  auto dataset = std::unique_ptr<H5::DataSet>{
    new H5::DataSet{file.openDataSet(name)}
  };

  // Verify that the dataset has the correct dimensionality and shape.
  std::array<hsize_t, 3> shape;
  auto ndim = dataset->getSpace().getSimpleExtentDims(shape.data());

  if (ndim != 3)
    throw std::invalid_argument{
      "dataset '" + name + "' in file '" + path + "' has " +
      std::to_string(ndim) + " dimensions (need 3)"
    };

  if (shape[2] != 3)
    throw std::invalid_argument{
      "dataset '" + name + "' in file '" + path + "' has " +
      std::to_string(shape[2]) + " columns (need 3)"
    };

  // Deduce number of configs and number of nucleons (A) from the shape.
  const auto& nconfigs = shape[0];
  const auto& A = shape[1];

  // Estimate the max radius from at least 500 nucleon positions.
  auto n = std::min(500/A + 1, nconfigs);
  std::array<hsize_t, 3> count = {n, A, 3};
  std::array<hsize_t, 3> start = {0, 0, 0};
  std::array<hsize_t, 2> shape_n = {n*A, 3};
  auto positions = read_dataset<float>(*dataset, count, start, shape_n);

  auto rmax_sq = 0.;

  for (const auto& position : positions) {
    auto& x = position[0];
    auto& y = position[1];
    auto& z = position[2];
    auto r_sq = x*x + y*y + z*z;
    if (r_sq > rmax_sq)
      rmax_sq = r_sq;
  }

  auto rmax = std::sqrt(rmax_sq);

  return std::unique_ptr<ManualNucleus>{
    new ManualNucleus{std::move(dataset), nconfigs, A, rmax}
  };
}

ManualNucleus::ManualNucleus(std::unique_ptr<H5::DataSet> dataset,
                             std::size_t nconfigs, std::size_t A, double rmax)
    : Nucleus(A),
      dataset_(std::move(dataset)),
      A_(A),
      rmax_(rmax),
      index_dist_(0, nconfigs - 1)
{}

ManualNucleus::~ManualNucleus() = default;

double ManualNucleus::radius() const {
  return rmax_;
}

void ManualNucleus::sample_nucleons_impl() {
  // Sample Euler rotation angles.
  // First is an azimuthal spin about the Z axis.
  const auto angle_1 = random::phi<double>();
  const auto c1 = std::cos(angle_1);
  const auto s1 = std::sin(angle_1);
  // Then a polar tilt about the original X axis, uniform in cos(theta).
  const auto c2 = random::cos_theta<double>();
  const auto s2 = std::sqrt(1. - c2*c2);
  // Finally another azimuthal spin about the original Z axis.
  const auto angle_3 = random::phi<double>();
  const auto c3 = std::cos(angle_3);
  const auto s3 = std::sin(angle_3);

  // Choose and read a random config from the dataset.
  std::array<hsize_t, 3> count = {1, A_, 3};
  std::array<hsize_t, 3> start = {index_dist_(random::engine), 0, 0};
  std::array<hsize_t, 2> shape = {A_, 3};
  const auto positions = read_dataset<float>(*dataset_, count, start, shape);

  // Loop over the positions and internal nucleons together.
  // Could use something fancy like boost::zip_iterator, but this is simpler and
  // more readable.
  auto positions_iter = positions.begin();
  for (auto&& nucleon : *this) {
    // Extract position vector and increment iterator.
    auto position = *positions_iter++;
    auto& x = position[0];
    auto& y = position[1];
    auto& z = position[2];

    // Rotate.
    auto x_rot = x*(c1*c3 - c2*s1*s3) - y*(c3*s1 + c1*c2*s3) + z*s2*s3;
    auto y_rot = x*(c1*s3 + c2*c3*s1) - y*(s1*s3 - c1*c2*c3) - z*c3*s2;

    set_nucleon_position(nucleon, x_rot, y_rot);
  }
}

#endif  // TRENTO_HDF5

}  // namespace trento
