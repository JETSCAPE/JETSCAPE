#ifndef LINE_H
#define LINE_H

#include <algorithm>
#include <array>

#include "GeneralGeometryElement.h"

namespace JetscapeCornelius {

/**
 * @class Line
 * @brief Represents a line in a geometric space, derived from
 * GeneralGeometryElement.
 *
 * The Line class encapsulates the properties and operations related to a line
 * segment in a geometric space. It provides methods to initialize the line,
 * flip its start and end points, and calculate various geometric properties
 * such as the normal and centroid.
 *
 * 01.10.2025 Hendrik Roch, Haydar Mehryar, Joe Latessa
 *
 */
class Line : public GeneralGeometryElement {
 protected:
  static constexpr int LINE_DIM =
      2;  ///< Dimension for line-specific properties
  static constexpr int LINE_CORNERS = 2;  ///< Number of corners for a line

  int x1, x2;       ///< Indices representing the line's dimensions
  int start_point;  ///< Index of the start point
  int end_point;    ///< Index of the end point
  std::array<std::array<double, DIM>, LINE_DIM>
      corners;                               ///< Array of line corners
  std::array<double, DIM> out;               ///< Output point of the line
  std::array<int, DIM - LINE_DIM> const_i;   ///< Constant indices for the line
  mutable std::array<double, DIM> reference_normal;  ///< Reference normal vector

 public:
  /**
   * @brief Default constructor for the Line class.
   *
   * Initializes a Line object. Calls the constructor of GeneralGeometryElement.
   */
  Line();

  /**
   * @brief Destructor for the Line class.
   *
   * Cleans up any resources used by the Line object.
   */
  ~Line();

  /**
   * @brief Initializes the line with specific corner points, an output point,
   * and constant indices.
   *
   * @param new_corners Array of corner points defining the line
   * @param new_out Point outside the surface
   * @param new_const_i Array of constant indices
   */
  inline void init_line(
      const std::array<std::array<double, DIM>, LINE_DIM>& new_corners,
      const std::array<double, DIM>& new_out,
      const std::array<int, DIM - LINE_DIM>& new_const_i) {
    corners = new_corners;
    out = new_out;
    const_i = new_const_i;
    start_point = 0;
    end_point = 1;
    // Use a lookup table for x1/x2 assignment
    constexpr int x_lookup[3][2][2] = {
        {{2, 3}, {1, 3}},  // new_const_i[0] == 0, new_const_i[1] == 1 or 2
        {{0, 3}, {0, 2}},  // new_const_i[0] == 1, new_const_i[1] == 2 or else
        {{0, 1}, {0, 1}}   // else
    };
    if (new_const_i[0] == 0) {
      x1 = (new_const_i[1] == 1) ? x_lookup[0][0][0] : x_lookup[0][1][0];
      x2 = (new_const_i[1] == 1) ? x_lookup[0][0][1] : x_lookup[0][1][1];
    } else if (new_const_i[0] == 1) {
      x1 = (new_const_i[1] == 2) ? x_lookup[1][0][0] : x_lookup[1][1][0];
      x2 = (new_const_i[1] == 2) ? x_lookup[1][0][1] : x_lookup[1][1][1];
    } else {
      x1 = x_lookup[2][0][0];
      x2 = x_lookup[2][0][1];
    }
    normal_calculated = centroid_calculated = false;
  }

  /**
   * @brief Flips the start and end points of the line.
   *
   * This method swaps the line's start and end points to reverse the direction
   * of the line.
   */
  inline void flip_start_end() { std::swap(start_point, end_point); }

  /**
   * @brief Calculates the normal vector of the line.
   *
   * Computes the normal vector for the line. This function must be implemented
   * based on the specific geometric context of the line.
   */
  inline void calculate_normal() const override {
    if (!centroid_calculated)
      calculate_centroid();
    // The normal is given by (-dy, dx)
    const double dx1 = corners[1][x1] - corners[0][x1];
    const double dx2 = corners[1][x2] - corners[0][x2];
    normal[x1] = -dx2;
    normal[x2] = dx1;
    normal[const_i[0]] = 0.0;
    normal[const_i[1]] = 0.0;
    for (int i = 0; i < DIM; ++i) {
      reference_normal[i] = out[i] - centroid[i];
    }
    flip_normal_if_needed(normal, reference_normal);
    normal_calculated = true;
  }

  /**
   * @brief Calculates the centroid of the line.
   *
   * Computes the centroid point of the line. This function must be implemented
   * based on the specific geometric context of the line.
   */
  inline void calculate_centroid() const override {
    for (int i = 0; i < DIM; i++) {
      centroid[i] = 0.5 * (corners[0][i] + corners[1][i]);
    }
    centroid_calculated = true;
  }

  /**
   * @brief Retrieves the start point of the line in 4D.
   *
   * @return Reference to the array representing the start point
   */
  inline const std::array<double, GeneralGeometryElement::DIM>&
  get_start_point() const {
    return corners[start_point];
  }

  /**
   * @brief Retrieves the end point of the line in 4D.
   *
   * @return Reference to the array representing the end point
   */
  inline const std::array<double, GeneralGeometryElement::DIM>& get_end_point()
      const {
    return corners[end_point];
  }

  /**
   * @brief Retrieves the point which is always outside in 4D.
   *
   * @return Reference to the array representing the outside point
   */
  inline const std::array<double, GeneralGeometryElement::DIM>&
  get_outside_point() const {
    return out;
  }
};

}  // namespace JetscapeCornelius

#endif  // LINE_H