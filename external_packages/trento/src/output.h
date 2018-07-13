// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// TRENTO3D: Three-dimensional extension of TRENTO by Weiyao Ke
// MIT License

#ifndef OUTPUT_H
#define OUTPUT_H

#include <functional>
#include <utility>
#include <vector>

#include "fwd_decl.h"

namespace trento {

/// Simple interface for outputting event data.  Determines which output formats
/// to create based on the configuration and writes those formats when called.
class Output {
 public:
  /// Instantiate from the configuration.
  Output(const VarMap& var_map);

  /// \rst
  /// Call the functor to output event data.  Arguments are perfect-forwarded to
  /// each output function.  The required arguments are
  ///
  /// - ``int`` event number
  /// - ``double`` impact parameter
  /// - ``const Event&`` Event object
  ///
  /// \endrst
  template <typename... Args>
  void operator()(Args&&... args) const;

 private:
  /// Internal storage of output functions.
  std::vector<std::function<void(int, double, const Event&)>> writers_;
};

template <typename... Args>
void Output::operator()(Args&&... args) const {
  for (const auto& write : writers_)
    write(std::forward<Args>(args)...);
}

}  // namespace trento

#endif  // OUTPUT_H
