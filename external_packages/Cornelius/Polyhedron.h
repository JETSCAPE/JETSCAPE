#ifndef POLYHEDRON_H
#define POLYHEDRON_H

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>

#include "GeneralGeometryElement.h"
#include "Line.h"
#include "Polygon.h"

namespace JetscapeCornelius {

/**
 * @class Polyhedron
 * @brief Represents a polyhedron composed of multiple polygons.
 *
 * The Polyhedron class inherits from GeneralGeometryElement and contains
 * methods for initializing the polyhedron, adding polygons, checking line
 * equality, calculating the volume of tetrahedrons, and calculating the
 * centroid and normal.
 *
 * 01.10.2025 Hendrik Roch, Haydar Mehryar, Joe Latessa
 *
 */
class Polyhedron : public GeneralGeometryElement {
 private:
  static constexpr int MAX_POLYGONS =
      24;  ///< Maximum number of polygons in a polyhedron
  static constexpr double INV_SIX = 1.0 / 6.0;  ///< Inverse of six.
  static constexpr double EPSILON = 1e-10;  ///< Epsilon value for comparison
    std::vector<const Polygon*>
      polygons;             ///< Pointers to polygons in the polyhedron
  int number_polygons;      ///< Number of polygons in the polyhedron
  int number_tetrahedrons;  ///< Number of tetrahedrons in the polyhedron

  // Arrays for temporary storage
  mutable std::array<double, DIM> Vout;  ///< Point outside the surface
  mutable std::array<double, DIM> a;     ///< Vector a of the triangle
  mutable std::array<double, DIM> b;     ///< Vector b of the triangle
  mutable std::array<double, DIM> c;     ///< Centroid of the triangle
  mutable std::array<double, DIM> n;     ///< Normal vector of the triangle
  mutable std::array<double, DIM> cm_i;  ///< Center of mass of the tetrahedron
  mutable std::array<std::array<double, DIM>, MAX_POLYGONS * 24>
      normals;  ///< Normal vectors

 public:
  /**
   * @brief Constructs a Polyhedron object.
   */
  Polyhedron();

  /**
   * @brief Destroys the Polyhedron object.
   */
  ~Polyhedron();

  /**
   * @brief Initializes the polyhedron.
   *
   * Sets the initial values for the indices, polygon count, tetrahedron count,
   * and flags.
   */
  inline void init_polyhedron() {
    // Reset the number of polygons and tetrahedrons in the polyhedron
    number_polygons = number_tetrahedrons = 0;
    // Set the flags for normal and centroid calculations to false
    normal_calculated = centroid_calculated = false;
    polygons.clear();
    polygons.reserve(MAX_POLYGONS);
    polygons.emplace_back(nullptr);  // Default placeholder
  }

  /**
   * @brief Adds a polygon to the polyhedron.
   *
   * @param new_polygon The polygon to add.
   * @param perform_no_check If true, adds the polygon without checking for
   * connection to existing polygons.
   * @return True if the polygon was added successfully, false otherwise.
   */
  bool add_polygon(const Polygon* new_polygon, bool perform_no_check);

  /**
   * @brief Checks if two lines are connected.
   *
   * @param line1 The first line.
   * @param line2 The second line.
   * @return True if the lines are connected, false otherwise.
   */
  inline bool lines_are_connected(const Line& line1, const Line& line2) {
    // Get the start and end points of the lines
    const auto& start_point1 = line1.get_start_point();
    const auto& end_point1 = line1.get_end_point();
    const auto& start_point2 = line2.get_start_point();

    // Check if any point pairs are close enough
    double difference1 = 0.0;
    double difference2 = 0.0;
    for (int i = 0; i < DIM; i++) {
      if (difference1 <= EPSILON)
        difference1 += std::abs(start_point1[i] - start_point2[i]);
      if (difference2 <= EPSILON)
        difference2 += std::abs(end_point1[i] - start_point2[i]);
      if (difference1 > EPSILON && difference2 > EPSILON) {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Calculates the normal vector of a tetrahedron, which is also
   * the volume of the tetrahedron.
   *
   * @param v1 First vector defining the tetrahedron.
   * @param v2 Second vector defining the tetrahedron.
   * @param v3 Third vector defining the tetrahedron.
   * @param n Normal vector of the tetrahedron.
   */
  inline void tetrahedron_volume(std::array<double, DIM>& v1,
                                 std::array<double, DIM>& v2,
                                 std::array<double, DIM>& v3,
                                 std::array<double, DIM>& n) const {
    // Calculate the volume of the tetrahedron
    const double bc01 = v2[0] * v3[1] - v2[1] * v3[0];
    const double bc02 = v2[0] * v3[2] - v2[2] * v3[0];
    const double bc03 = v2[0] * v3[3] - v2[3] * v3[0];
    const double bc12 = v2[1] * v3[2] - v2[2] * v3[1];
    const double bc13 = v2[1] * v3[3] - v2[3] * v3[1];
    const double bc23 = v2[2] * v3[3] - v2[3] * v3[2];
    n[0] = (v1[1] * bc23 - v1[2] * bc13 + v1[3] * bc12) * INV_SIX;
    n[1] = -(v1[0] * bc23 - v1[2] * bc03 + v1[3] * bc02) * INV_SIX;
    n[2] = (v1[0] * bc13 - v1[1] * bc03 + v1[3] * bc01) * INV_SIX;
    n[3] = -(v1[0] * bc12 - v1[1] * bc02 + v1[2] * bc01) * INV_SIX;
  }

  /**
   * @brief Calculates the centroid of the polyhedron.
   *
   * Computes the centroid as the volume-weighted average of the individual
   * tetrahedrons.
   */
  void calculate_centroid() const override;

  /**
   * @brief Calculates the normal of the polyhedron.
   *
   * Computes the normal as the sum of the normals of the individual
   * tetrahedrons.
   */
  void calculate_normal() const override;

  /**
   * @brief Retrieves the number of polygons in the polyhedron.
   *
   * @return The number of polygons in the polyhedron.
   */
  inline int get_number_polygons() const { return number_polygons; }

  /**
   * @brief Retrieves the number of tetrahedrons in the polyhedron.
   *
   * @return The number of tetrahedrons in the polyhedron.
   */
  inline int get_number_tetrahedrons() const { return number_tetrahedrons; }
};

}  // namespace JetscapeCornelius

#endif  // POLYHEDRON_H