#ifndef POLYGON_H
#define POLYGON_H

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <numeric>
#include <vector>

#include "GeneralGeometryElement.h"
#include "Line.h"

namespace JetscapeCornelius {

/**
 * @class Polygon
 * @brief A class representing a polygon in a geometric space.
 *
 * This class extends the GeneralGeometryElement class and provides methods
 * to manage and calculate properties related to polygons.
 *
 * 01.10.2025 Hendrik Roch, Haydar Mehryar, Joe Latessa
 *
 */
class Polygon : public GeneralGeometryElement {
 protected:
    static constexpr int MAX_LINES =
      24;  ///< Maximum number of lines in a polygon
    std::vector<const Line*> lines;  ///< Pointers to lines in the polygon
    int number_lines;                 ///< Number of lines in the polygon
  int x1, x2, x3;           ///< Indices representing the polygon's dimensions
  int const_i;              ///< Constant index for the polygon

  // Arrays a and b to store the vectors of the triangles
  mutable std::array<double, DIM> a;                  ///< Vector a of the triangle
  mutable std::array<double, DIM> b;                  ///< Vector b of the triangle
  mutable std::array<double, DIM> triangle_centroid;  ///< Centroid of the triangle

  mutable std::array<std::array<double, DIM>, MAX_LINES> normals;  ///< Normal vectors

  static constexpr double EPSILON = 1e-10;  ///< Small value for epsilon.

 public:
  /**
   * @brief Default constructor for the Polygon class.
   *
   * Initializes a new instance of the Polygon class.
   */
  Polygon();

  /**
   * @brief Destructor for the Polygon class.
   *
   * Cleans up resources used by the Polygon class.
   */
  ~Polygon();

  /**
   * @brief Initializes the polygon with a constant index.
   *
   * @param new_const_i The constant index for the polygon's dimensions.
   */
  inline void init_polygon(int new_const_i) {
    // Copy the new value into the class variable
    const_i = new_const_i;
    // Fix the indices which are not constant
    switch (new_const_i) {
      case 0:
        x1 = 1;
        x2 = 2;
        x3 = 3;
        break;
      case 1:
        x1 = 0;
        x2 = 2;
        x3 = 3;
        break;
      case 2:
        x1 = 0;
        x2 = 1;
        x3 = 3;
        break;
      case 3:
        x1 = 0;
        x2 = 1;
        x3 = 2;
        break;
      default:
        break;
    }
    // Set the flags for normal and centroid calculations to false
    normal_calculated = centroid_calculated = false;
    // Reset the number of lines in the polygon
    number_lines = 0;
    lines.clear();
    lines.reserve(MAX_LINES);
  }

  /**
   * @brief Adds a line to the polygon.
   *
   * @param new_line The line to be added.
   * @param perform_no_check If true, the line is added without connectivity
   * checks.
   * @return True if the line was successfully added, otherwise false.
   */
  bool add_line(const Line* new_line, bool perform_no_check);

  /**
   * @brief Gets the number of lines in the polygon.
   *
   * @return The number of lines in the polygon.
   */
  inline int get_number_lines() const { return number_lines; }

  /**
   * @brief Calculates the normal vector of the polygon.
   *
   * This method calculates the normal vector based on the lines and centroid of
   * the polygon.
   */
  void calculate_normal() const override;

  /**
   * @brief Calculates the centroid of the polygon.
   *
   * This method calculates the centroid based on the vertices of the polygon.
   */
  void calculate_centroid() const override;

  /**
   * @brief Gets the lines that form the polygon.
   *
   * @return A reference to the array of lines in the polygon.
   */
  inline const std::vector<const Line*>& get_lines() const { return lines; }

  /**
   * @brief Prints the triangles formed from the polygon into a given file.
   * Prints the absolute points, so this file can be used to plot the surface.
   *
   * @param file The output file stream.
   * @param position The position offset for printing the polygon.
   */
  void print(std::ofstream& file, std::array<double, DIM> position);
};

}  // namespace JetscapeCornelius

#endif  // POLYGON_H