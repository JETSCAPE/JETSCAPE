// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// TRENTO3D: Three-dimensional extension of TRENTO by Weiyao Ke
// MIT License

#ifndef NUCLEUS_H
#define NUCLEUS_H

#include <memory>
#include <random>
#include <string>
#include <vector>

#include "fwd_decl.h"
#include "nucleon.h"

#ifdef TRENTO_HDF5
// forward declaration for std::unique_ptr<H5::DataSet> member of ManualNucleus
namespace H5 {
class DataSet;
}
#endif

namespace trento {

// Alias for a smart pointer to a Nucleus.
using NucleusPtr = std::unique_ptr<Nucleus>;

/// \rst
/// Interface class to all nucleus types.  Stores an ensemble of nucleons and
/// randomly samples their positions.  Implements a standard iterator interface
/// through ``begin()`` and ``end()`` functions.  Iterating over a ``Nucleus``
/// means iterating over its ``Nucleon`` members.
/// \endrst
class Nucleus {
 public:
  /// \rst
  /// The canonical way to create a ``Nucleus``.  Constructs the appropriate
  /// subclass for the given species.  Sets the Woods-Saxon parameters for Au,
  /// Pb, U, etc; parameters copied from the PHOBOS Glauber model:
  ///
  /// - http://inspirehep.net/record/786828
  /// - http://inspirehep.net/record/1310629
  ///
  /// Example::
  ///
  ///   std::unique_ptr<Nucleus> lead_nucleus = Nucleus::create("Pb");
  ///   double radius = lead_nucleus->radius();
  ///   lead_nucleus->sample_nucleons(0);
  ///   for (const auto& nucleon : *lead_nucleus)
  ///     do_something(nucleon)
  ///
  /// \endrst
  ///
  /// \param species standard symbol, e.g. "p" for proton or "Pb" for lead-208
  /// \param nucleon_dmin minimum nucleon-nucleon distance for Woods-Saxon
  /// nuclei (optional, default zero)
  ///
  /// \return a smart pointer \c std::unique_ptr<Nucleus>
  ///
  /// \throw std::invalid_argument for unknown species
  static NucleusPtr create(const std::string& species, double nucleon_width, double nucleon_dmin = 0);

  /// Default virtual destructor for abstract base class.
  virtual ~Nucleus() = default;

  /// The "radius", i.e. the maximum distance at which a nucleon could be
  /// placed.
  virtual double radius() const = 0;

  /// Sample a new ensemble of nucleon positions with the given offset in the
  /// x-direction.
  void sample_nucleons(double offset);

  using size_type = std::vector<Nucleon>::size_type;
  using iterator = std::vector<Nucleon>::iterator;
  using const_iterator = std::vector<Nucleon>::const_iterator;

  // size = number of nucleons
  size_type size() const noexcept
  { return nucleons_.size(); }

  // non-const overload
  iterator begin() noexcept
  { return nucleons_.begin(); }
  iterator end() noexcept
  { return nucleons_.end(); }

  // const overload
  const_iterator begin() const noexcept
  { return nucleons_.begin(); }
  const_iterator end() const noexcept
  { return nucleons_.end(); }

  // forced const
  const_iterator cbegin() const noexcept
  { return nucleons_.cbegin(); }
  const_iterator cend() const noexcept
  { return nucleons_.cend(); }

 protected:
  /// Constructor only accessible by derived classes.
  /// \param A number of nucleons
  explicit Nucleus(std::size_t A);

  /// \rst
  /// Set a ``Nucleon`` position.  This function provides a consistent interface
  /// to derived classes and ensures the position is correctly offset.
  /// ``Nucleus`` is a friend of ``Nucleon`` and therefore able to set nucleon
  /// positions; the derived classes must use this function to set positions.
  /// \endrst
  void set_nucleon_position(iterator nucleon, double x, double y, double z);

 private:
  /// Internal interface to the actual implementation of the nucleon sampling
  /// algorithm, used in public function sample_nucleons().  This function must
  /// sample nucleon positions relative to the origin and set them using the
  /// protected function set_nucleon_position(), which enforces the offset.
  virtual void sample_nucleons_impl() = 0;

  /// Internal storage of Nucleon objects.
  std::vector<Nucleon> nucleons_;

  /// Offset of nucleon x-positions.
  /// This variable is reset upon each call of sample_nucleons() and is read by
  /// set_nucleon_position().
  double offset_;
};

// Now declare Nucleus subclasses.
// Each subclass must do the following:
//   - Define a public constructor which (at minimum) initializes a Nucleus
//     with the required number of nucleons.
//   - Override and implement the pure virtual functions
//     radius() and sample_nucleons_impl().

/// Trivial nucleus with a single nucleon.
class Proton : public Nucleus {
 public:
  /// Default constructor.
  Proton();

  /// The radius of a proton is trivially zero.
  virtual double radius() const override;

 private:
  /// The proton trivially places its nucleon at the origin.
  virtual void sample_nucleons_impl() override;
};

/// \rst
/// Samples pairs of nucleons from the Hulthén wavefunction
///
/// .. math::
///
///   f(r) \propto \biggl( \frac{\exp(-ar) - \exp(-br)}{r} \biggr)^2.
///
/// http://inspirehep.net/record/1261055
///
/// \endrst
class Deuteron : public Nucleus {
 public:
  /// Insantiate with values for the (a, b) parameters.
  /// The defaults are from the PHOBOS Glauber model.
  Deuteron(double a = 0.457, double b = 2.35);

  /// The radius is computed from the parameters (a, b).
  virtual double radius() const override;

 private:
  /// Sample positions from the Hulthén wavefunction.
  virtual void sample_nucleons_impl() override;

  /// Internal storage of wavefunction parameters.
  const double a_, b_;
};

/// A nucleus that can check if its nucleons are within a specified minimum
/// distance.
class MinDistNucleus : public Nucleus {
 protected:
  /// \param A number of nucleons
  /// \param dmin minimum nucleon-nucleon distance (optional, default zero)
  MinDistNucleus(std::size_t A, double dmin = 0);

  /// \rst
  /// Check if a ``Nucleon`` is too close (within the minimum distance) of any
  /// previously placed nucleons.  Specifically, check nucleons from ``begin()``
  /// up to the given iterator.
  /// \endrst
  bool is_too_close(const_iterator nucleon) const;

 private:
  /// Internal storage of squared minimum distance.
  const double dminsq_;
};

/// \rst
/// Samples nucleons from a spherically symmetric Woods-Saxon distribution
///
/// .. math::
///
///   f(r) \propto \frac{1}{1 + \exp(\frac{r-R}{a})}.
///
/// For non-deformed heavy nuclei such as lead.
///
/// \endrst
class WoodsSaxonNucleus : public MinDistNucleus {
 public:
  /// ``Nucleus::create()`` sets these parameters for a given species.
  /// \param A number of nucleons
  /// \param R Woods-Saxon radius
  /// \param a Woods-Saxon surface thickness
  /// \param dmin minimum nucleon-nucleon distance (optional, default zero)
  WoodsSaxonNucleus(std::size_t A, double R, double a, double dmin = 0);

  /// The radius of a Woods-Saxon Nucleus is computed from the parameters (R, a).
  virtual double radius() const override;

 private:
  /// Sample Woods-Saxon nucleon positions.
  virtual void sample_nucleons_impl() override;

  /// Woods-Saxon parameters.
  const double R_, a_;

  /// Woods-Saxon distribution object.  Since the dist does not have an analytic
  /// inverse CDF, approximate it as a piecewise linear dist.  For a large
  /// number of steps this is very accurate.
  mutable std::piecewise_linear_distribution<double> woods_saxon_dist_;
};

/// \rst
/// Samples nucleons from a deformed spheroidal Woods-Saxon distribution
///
/// .. math::
///
///   f(r, \theta) \propto
///   \frac{1}{1 + \exp(\frac{r-R(1+\beta_2Y_{20}+\beta_4Y_{40})}{a})}.
///
/// For deformed heavy nuclei such as uranium.
///
/// \endrst
class DeformedWoodsSaxonNucleus : public MinDistNucleus {
 public:
  /// ``Nucleus::create()`` sets these parameters for a given species.
  /// \param A number of nucleons
  /// \param R Woods-Saxon radius
  /// \param a Woods-Saxon surface thickness
  /// \param beta2 Woods-Saxon deformation parameter
  /// \param beta4 Woods-Saxon deformation parameter
  /// \param dmin minimum nucleon-nucleon distance (optional, default zero)
  DeformedWoodsSaxonNucleus(std::size_t A, double R, double a,
                            double beta2, double beta4, double dmin = 0);

  /// The radius of a deformed Woods-Saxon Nucleus is computed from the
  /// parameters (R, a, beta2, beta4).
  virtual double radius() const override;

 private:
  /// Sample deformed Woods-Saxon nucleon positions.
  virtual void sample_nucleons_impl() override;

  /// Evaluate the deformed Woods-Saxon distribution.
  double deformed_woods_saxon_dist(double r, double cos_theta) const;

  /// Woods-Saxon parameters.
  const double R_, a_, beta2_, beta4_;

  /// Maximum radius.
  const double rmax_;
};

#ifdef TRENTO_HDF5

/// Reads manual nuclear configurations from an HDF5 file.
class ManualNucleus : public Nucleus {
 public:
  /// Create a ManualNucleus that reads from the given file.
  /// Throw std::invalid_argument if there are any problems.
  ///
  /// Since this is a derived class, the base Nucleus class is initialized
  /// before any the data members.  This creates a bit of a catch-22, since the
  /// number of nucleons must be known to initialize the base class, but that
  /// must be deduced from the file.  As a workaround, this factory function
  /// opens the file, determines the number of nucleons, and then calls the
  /// constructor.
  static std::unique_ptr<ManualNucleus> create(const std::string& path);

  /// Must define destructor because of member pointer to incomplete type.
  /// See explanation for Collider destructor.
  virtual ~ManualNucleus() override;

  /// The radius is determined by reading many positions from the file and
  /// saving the maximum.
  virtual double radius() const override;

 private:
  /// Private constructor -- use create().
  /// \param dataset smart pointer to HDF5 dataset
  /// \param nconfigs number of nucleus configs in the dataset
  /// \param A number of nucleons
  /// \param rmax max radius
  ManualNucleus(std::unique_ptr<H5::DataSet> dataset,
                std::size_t nconfigs, std::size_t A, double rmax);

  /// Read a configuration from the file, rotate it, and set nucleon positions.
  virtual void sample_nucleons_impl() override;

  /// Internal pointer to HDF5 dataset object (PIMPL-like).
  const std::unique_ptr<H5::DataSet> dataset_;

  /// Internal storage of the maximum radius.
  const double rmax_;

  /// Distribution for choosing random configs.
  std::uniform_int_distribution<std::size_t> index_dist_;
};

#endif  // TRENTO_HDF5

}  // namespace trento

#endif  // NUCLEUS_H
