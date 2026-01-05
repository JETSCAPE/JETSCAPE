#include "Hypercube.h"

namespace JetscapeCornelius {

Hypercube::Hypercube() : number_polyhedra(0), ambiguous(false) {
  polyhedra.clear();
  polyhedra.reserve(MAX_POLYHEDRONS);
}

Hypercube::~Hypercube() = default;

int Hypercube::split_to_cubes(double value) {
  int number_points_below_value = 0;
  int cube_index = 0;
  // Use std::array for c_v
  for (int i = 0; i < DIM; ++i) {
    const int c_i = i;
    for (int j = 0; j < STEPS; ++j) {
      std::array<double, DIM> c_v = {0.0, 0.0, 0.0, 0.0};
      c_v[i] = j * dx[i];
      for (int ci1 = 0; ci1 < STEPS; ++ci1) {
        for (int ci2 = 0; ci2 < STEPS; ++ci2) {
          for (int ci3 = 0; ci3 < STEPS; ++ci3) {
            double hypercube_value;
            if (i == 0) {
              hypercube_value = hypercube[j][ci1][ci2][ci3];
              if (hypercube_value < value) {
                ++number_points_below_value;
              }
            } else if (i == 1) {
              hypercube_value = hypercube[ci1][j][ci2][ci3];
            } else if (i == 2) {
              hypercube_value = hypercube[ci1][ci2][j][ci3];
            } else {
              hypercube_value = hypercube[ci1][ci2][ci3][j];
            }
            cube[ci1][ci2][ci3] = hypercube_value;
          }
        }
      }
      cubes[cube_index++].init_cube(cube, c_i, c_v[i], dx);
    }
  }
  return number_points_below_value;
}

void Hypercube::construct_polyhedra(double value) {
  const int number_points_below_value = split_to_cubes(value);

  // Store the reference to the polygons
  int number_polygons = 0;
  for (int i = 0; i < NCUBES; i++) {
    cubes[i].construct_polygons(value);
    auto& polygons_cube = cubes[i].get_polygons();
    const int n_local = cubes[i].get_number_polygons();
    for (int j = 0; j < n_local; ++j) {
      polygon_refs[number_polygons++] = &polygons_cube[j];
    }
  }

  check_ambiguity(number_points_below_value);

  if (number_polygons == 0) {
    return;
  }

  if (ambiguous) {
    // The surface might be ambiguous and we need to connect the polygons and
    // see how many polyhedra we have
    not_used.fill(1);
    // Keep track of the used number of lines
    int used = 0;
    do {
      // Ensure there's space in the vector
      if (number_polyhedra >= polyhedra.size()) {
        polyhedra.emplace_back();  // Add a new Polyhedron if needed
      }
      polyhedra[number_polyhedra].init_polyhedron();
      // Go through all the polygons and try to add them to the polyhedron
      for (int i = 0; i < number_polygons; i++) {
        if (not_used[i]) {
          if (polyhedra[number_polyhedra].add_polygon(polygon_refs[i],
                                                      false)) {
            not_used[i] = 0;
            used++;
            // If the polygon is successfully added we start the loop from the
            // beginning
            i = 0;
          }
        }
      }
      number_polyhedra++;
    } while (used < number_polygons);
  } else {
    // Here surface cannot be ambiguous and all polygons can be added to
    // the polyhedron without ordering them
    if (number_polyhedra >= polyhedra.size()) {
      polyhedra.emplace_back();  // Add a new Polyhedron if needed
    }
    polyhedra[number_polyhedra].init_polyhedron();
    for (int i = 0; i < number_polygons; i++) {
      polyhedra[number_polyhedra].add_polygon(polygon_refs[i], true);
    }
    number_polyhedra++;
  }
}

}  // namespace JetscapeCornelius
