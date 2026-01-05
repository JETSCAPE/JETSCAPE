#include "Cube.h"

namespace JetscapeCornelius {

Cube::Cube() : number_lines(0), number_polygons(0), ambiguous(false) {
  polygons.clear();
  polygons.reserve(MAX_POLYGONS);
  polygons.emplace_back();  // Default to construct 1 Polygon
}

Cube::~Cube() = default;

void Cube::split_to_squares() {
  int number_squares = 0;
  for (int i = 0; i < DIM; ++i) {
    if (i == const_i)
      continue;

    std::array<int, STEPS> c_i = {const_i, i};
    for (int j = 0; j < STEPS; ++j) {
      std::array<double, STEPS> c_v = {const_value, j * dx[i]};

      // Optimize all cases with memcpy if possible
      if (i == x1) {
        for (int ci1 = 0; ci1 < STEPS; ++ci1) {
          std::memcpy(&square[ci1][0], &cube[j][ci1][0],
                      STEPS * sizeof(double));
        }
      } else if (i == x2) {
        for (int ci1 = 0; ci1 < STEPS; ++ci1) {
          std::memcpy(&square[ci1][0], &cube[ci1][j][0],
                      STEPS * sizeof(double));
        }
      } else {  // i == x3
        for (int ci1 = 0; ci1 < STEPS; ++ci1) {
          for (int ci2 = 0; ci2 < STEPS; ++ci2) {
            square[ci1][ci2] = cube[ci1][ci2][j];
          }
        }
      }
      squares[number_squares++].init_square(square, c_i, c_v, dx);
    }
  }
}

void Cube::construct_polygons(double value) {
  // Start by splitting the cube to squares and finding the lines
  split_to_squares();

  // Store the reference to the lines
  number_lines = 0;
  for (int i = 0; i < NSQUARES; i++) {
    squares[i].construct_lines(value);
    auto& lines_square = squares[i].get_lines();
    const int n_local = squares[i].get_number_lines();
    for (int j = 0; j < n_local; ++j) {
      line_refs[number_lines++] = &lines_square[j];
    }
  }

  // If no lines were found we may exit. This can happen only in 4D case
  if (number_lines == 0) {
    return;
  }
  // Then we check if the surface is ambiguous and continue
  check_ambiguity(number_lines);
  if (ambiguous) {
    // Surface is ambiguous, connect the lines to polygons and see how
    // many polygons we have
    not_used.fill(1);
    // Keep track of the lines which are used
    int used = 0;
    do {
      if (number_lines - used < 3) {
        std::cerr << "Error: cannot construct polygon from "
                  << number_lines - used << " lines" << std::endl;
        exit(1);
      }
      // Ensure there's space in the vector
      if (number_polygons >= polygons.size()) {
        polygons.emplace_back();  // Add a new Polygon if needed
      }
      // Initialize a new polygon
      polygons[number_polygons].init_polygon(const_i);
      // Go through all lines and try to add them to the polygon
      for (int i = 0; i < number_lines; i++) {
        // add_line returns true if line is successfully added
        if (not_used[i]) {
          if (polygons[number_polygons].add_line(line_refs[i], false)) {
            not_used[i] = 0;
            used++;
            // If line is successfully added we start the loop from the
            // beginning
            i = 0;
          }
        }
      }
      // When we have reached this point one complete polygon is formed
      number_polygons++;
    } while (used < number_lines);
  } else {
    // Surface is not ambiguous, so we have only one polygon and all lines
    // can be added to it without ordering them
    if (number_polygons >= polygons.size()) {
      polygons.emplace_back();  // Add a new Polygon if needed
    }
    polygons[number_polygons].init_polygon(const_i);
    for (int i = 0; i < number_lines; i++) {
      polygons[number_polygons].add_line(line_refs[i], true);
    }
    number_polygons++;
  }
}

}  // namespace JetscapeCornelius