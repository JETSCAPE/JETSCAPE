#ifndef GENERAL_GEOMETRY_ELEMENT_H
#define GENERAL_GEOMETRY_ELEMENT_H

#include <algorithm>
#include <array>
#include <numeric>

namespace JetscapeCornelius {

/**
 * @class GeneralGeometryElement
 * @brief A class representing a general geometric element with normal and
 * centroid calculations.
 *
 * This class provides functionality to manage and manipulate geometric
 * elements, including methods for calculating normals and centroids, and
 * flipping the orientation of the normal if needed.
 *
 * 01.10.2025 Hendrik Roch, Haydar Mehryar, Joe Latessa
 *
 */
class GeneralGeometryElement {
 protected:
  /// Dimension of the geometric space (e.g., 4 for 4D space).
  static constexpr int DIM = 4;

  /// Flag indicating if the normal has been calculated.
  mutable bool normal_calculated;

  /// Flag indicating if the centroid has been calculated.
  mutable bool centroid_calculated;

  /// Array representing the normal vector of the geometric element.
  mutable std::array<double, DIM> normal;

  /// Array representing the centroid of the geometric element.
  mutable std::array<double, DIM> centroid;

 public:
  /**
   * @brief Constructs a GeneralGeometryElement object.
   *
   * Initializes the flags for normal and centroid calculations to false.
   */
  GeneralGeometryElement();

  /**
   * @brief Destroys the GeneralGeometryElement object.
   *
   * Virtual destructor to allow proper cleanup in derived classes.
   */
  ~GeneralGeometryElement();

  /**
   * @brief Calculates the normal vector of the geometric element.
   *
   * This is a virtual method meant to be overridden by derived classes.
   * The default implementation does nothing.
   */
  virtual void calculate_normal() const;

  /**
   * @brief Calculates the centroid of the geometric element.
   *
   * This is a virtual method meant to be overridden by derived classes.
   * The default implementation does nothing.
   */
  virtual void calculate_centroid() const;

  /**
   * @brief Gets the normal vector of the geometric element.
   *
   * @return A constant array representing the normal vector.
   */
  inline std::array<double, DIM>& get_normal() const {
    if (!normal_calculated)
      calculate_normal();
    return normal;
  }

  /**
   * @brief Gets the centroid of the geometric element.
   *
   * @return A constant array representing the centroid.
   */
  inline std::array<double, DIM>& get_centroid() const {
    if (!centroid_calculated)
      calculate_centroid();
    return centroid;
  }

  /**
   * @brief Adjusts the normal vector if needed based on a reference normal.
   *
   * If the dot product between the provided normal and the reference normal is
   * negative, the normal vector is flipped.
   *
   * @param normal The normal vector to be checked and maybe adjusted.
   * @param reference_normal The reference normal vector used for comparison.
   */
  inline void flip_normal_if_needed(std::array<double, DIM>& normal,
                                    std::array<double, DIM>& reference_normal)
      const {
    const double dot_product = std::inner_product(
        normal.begin(), normal.end(), reference_normal.begin(), 0.0);
    if (dot_product < 0) {
      std::transform(normal.begin(), normal.end(), normal.begin(),
                     std::negate<double>());
    }
  }
};

}  // namespace JetscapeCornelius

#endif  // GENERAL_GEOMETRY_ELEMENT_H