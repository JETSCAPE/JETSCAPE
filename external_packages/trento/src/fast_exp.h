// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// TRENTO3D: Three-dimensional extension of TRENTO by Weiyao Ke
// MIT License

#ifndef FAST_EXP_H
#define FAST_EXP_H

#include <cmath>
#include <stdexcept>
#include <vector>

namespace trento {

/// \rst
/// Fast exponential approximation, to be used as a drop-in replacement for
/// ``std::exp`` when it will be evaluated many times within a fixed range.
/// Works by pre-tabulating exp() values and exploiting the leading-order Taylor
/// expansion; for step size `dx` the error is `\mathcal O(dx^2)`.
///
/// Example::
///
///   FastExp<double> fast_exp{0., 1., 11};  // tabulate at 0.0, 0.1, 0.2, ...
///   fast_exp(0.50);  // evaluate at table point -> exact result
///   fast_exp(0.55);  // midway between points -> worst-case error
///
/// \endrst
template <typename T = double>
class FastExp {
 public:
  /// Pre-tabulate exp() values from \c xmin to \c xmax in \c nsteps
  /// evenly-spaced intervals.
  FastExp(T xmin, T xmax, std::size_t nsteps);

  /// Evaluate the exponential at \em x (must be within range).
  T operator()(T x) const;

 private:
  /// Minimum and maximum.
  const T xmin_, xmax_;

  /// Step size.
  const T dx_;

  /// Tabulated exp() values.
  std::vector<T> table_;
};

template <typename T>
FastExp<T>::FastExp(T xmin, T xmax, std::size_t nsteps)
    : xmin_(xmin),
      xmax_(xmax),
      dx_((xmax-xmin)/(nsteps-1)),
      table_(nsteps) {
  // Tabulate evenly-spaced exp() values.
  for (std::size_t i = 0; i < nsteps; ++i)
    table_[i] = std::exp(xmin_ + i*dx_);
}

template <typename T>
inline T FastExp<T>::operator()(T x) const {
#ifndef NDEBUG
  if (x < xmin_ || x > xmax_)
    throw std::out_of_range{"argument must be within [xmin, xmax]"};
#endif

  // Determine the table index of the nearest tabulated value.
  auto index = static_cast<std::size_t>((x - xmin_)/dx_ + .5);

  // Compute the leading-order Taylor expansion.
  // exp(x) = exp(x0) * exp(x-x0) =~ exp(x0) * (1 + x - x0)
  // exp(x0) = table_[index]
  // x0 = xmin_ + index*dx_
  return table_[index] * (1. + x - xmin_ - index*dx_);
}

}  // namespace trento

#endif  // FAST_EXP_H
