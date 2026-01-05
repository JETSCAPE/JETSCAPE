#include "Polyhedron.h"

#include <iostream>

namespace JetscapeCornelius {

Polyhedron::Polyhedron() {
  polygons.clear();
  polygons.reserve(MAX_POLYGONS);
  polygons.emplace_back(nullptr);  // Default placeholder
}

Polyhedron::~Polyhedron() = default;

bool Polyhedron::add_polygon(const Polygon* new_polygon,
                             bool perform_no_check) {
  // For the first polygon, we don't need to check
  if (number_polygons == 0 || perform_no_check) {
    if (number_polygons >= polygons.size()) {
      polygons.emplace_back(nullptr);
    }
                      polygons[number_polygons++] = new_polygon;  // Assign new polygon
    number_tetrahedrons += new_polygon->get_number_lines();
    return true;
  }
  const int number_lines1 = new_polygon->get_number_lines();
  const auto& new_lines = new_polygon->get_lines();
  for (int i = 0; i < number_polygons; ++i) {
    const int n_lines2 = polygons[i]->get_number_lines();
    const auto& lines2 = polygons[i]->get_lines();
    for (int j = 0; j < number_lines1; ++j) {
      for (int k = 0; k < n_lines2; ++k) {
        if (lines_are_connected(*new_lines[j], *lines2[k])) {
          if (number_polygons >= polygons.size()) {
            polygons.emplace_back(nullptr);
          }
          polygons[number_polygons++] = new_polygon;
          number_tetrahedrons += number_lines1;
          return true;
        }
      }
    }
  }
  return false;
}

void Polyhedron::calculate_centroid() const {
  // Array of 0s to store the mean values
  std::array<double, DIM> mean_values = {0};
  for (int i = 0; i < number_polygons; ++i) {
    const auto& lines = polygons[i]->get_lines();
    int n_lines = polygons[i]->get_number_lines();
    for (int j = 0; j < n_lines; ++j) {
      const auto& start_point = lines[j]->get_start_point();
      const auto& end_point = lines[j]->get_end_point();
      for (int k = 0; k < DIM; ++k) {
        mean_values[k] += start_point[k] + end_point[k];
      }
    }
  }
  for (int k = 0; k < DIM; ++k) {
    mean_values[k] /= (2.0 * number_tetrahedrons);
  }
  std::array<double, DIM> sum_up = {0};
  double sum_down = 0.0;
  for (int i = 0; i < number_polygons; ++i) {
    const auto& lines = polygons[i]->get_lines();
    int n_lines = polygons[i]->get_number_lines();
    const auto& cent = polygons[i]->get_centroid();
    for (int j = 0; j < n_lines; ++j) {
      const auto& start_point = lines[j]->get_start_point();
      const auto& end_point = lines[j]->get_end_point();
      for (int k = 0; k < DIM; ++k) {
        cm_i[k] =
            (start_point[k] + end_point[k] + cent[k] + mean_values[k]) * 0.25;
        a[k] = start_point[k] - mean_values[k];
        b[k] = end_point[k] - mean_values[k];
        c[k] = cent[k] - mean_values[k];
      }
      tetrahedron_volume(a, b, c, n);
      double V_i = 0.0;
      for (int k = 0; k < DIM; ++k) {
        V_i += n[k] * n[k];
      }
      V_i = std::sqrt(V_i);
      for (int k = 0; k < DIM; ++k) {
        sum_up[k] += V_i * cm_i[k];
      }
      sum_down += V_i;
    }
  }
  for (int k = 0; k < DIM; ++k) {
    centroid[k] = sum_up[k] / sum_down;
  }
  centroid_calculated = true;
}

void Polyhedron::calculate_normal() const {
  if (!centroid_calculated) {
    calculate_centroid();
  }
  int index_tetrahedron = 0;
  for (int i = 0; i < number_polygons; ++i) {
    const auto& lines = polygons[i]->get_lines();
    const auto& cent = polygons[i]->get_centroid();
    int n_lines = polygons[i]->get_number_lines();
    for (int j = 0; j < n_lines; ++j) {
      const auto& start_point = lines[j]->get_start_point();
      const auto& end_point = lines[j]->get_end_point();
      const auto& o = lines[j]->get_outside_point();
      for (int k = 0; k < DIM; ++k) {
        a[k] = start_point[k] - centroid[k];
        b[k] = end_point[k] - centroid[k];
        c[k] = cent[k] - centroid[k];
        Vout[k] = o[k] - centroid[k];
      }
      tetrahedron_volume(a, b, c, normals[index_tetrahedron]);
      flip_normal_if_needed(normals[index_tetrahedron], Vout);
      ++index_tetrahedron;
    }
  }
  for (int k = 0; k < DIM; ++k) {
    normal[k] = 0.0;
    for (int j = 0; j < number_tetrahedrons; ++j) {
      normal[k] += normals[j][k];
    }
  }
  normal_calculated = true;
}

}  // namespace JetscapeCornelius