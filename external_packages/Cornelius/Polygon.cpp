#include "Polygon.h"

#include <iostream>

namespace JetscapeCornelius {

Polygon::Polygon() {
  lines.clear();
  lines.reserve(MAX_LINES);
}

Polygon::~Polygon() = default;

bool Polygon::add_line(const Line* new_line, bool perform_no_check) {
  // std::cout << "Polygon::add_line, size = " << lines.size() << ", " <<
  // number_lines << std::endl;
  //  For the first line, we don't need to check
  if (number_lines == 0 || perform_no_check) {
    if (number_lines >= lines.size()) {
      lines.emplace_back();
    }
    lines[number_lines++] = new_line;
    return true;
  }
  const auto& start_point = new_line->get_start_point();
  const auto& end_point = new_line->get_end_point();
  const auto& last_end_point = lines[number_lines - 1]->get_end_point();
  double difference1 = 0.0, difference2 = 0.0;
  for (int i = 0; i < DIM; ++i) {
    difference1 += std::abs(start_point[i] - last_end_point[i]);
    difference2 += std::abs(end_point[i] - last_end_point[i]);
  }
  if (difference1 < EPSILON || difference2 < EPSILON) {
    if (difference2 < EPSILON) {
      const_cast<Line*>(new_line)->flip_start_end();
    }
    if (number_lines >= lines.size()) {
      lines.emplace_back();
    }
    lines[number_lines++] = new_line;
    return true;
  }
  return false;
}

void Polygon::calculate_centroid() const {
  // Array of 0s to store the mean values
  std::array<double, DIM> mean_values = {0};
  for (int i = 0; i < number_lines; ++i) {
    const auto& start_point = lines[i]->get_start_point();
    const auto& end_point = lines[i]->get_end_point();
    for (int k = 0; k < DIM; ++k) {
      mean_values[k] += start_point[k] + end_point[k];
    }
  }
  for (int i = 0; i < DIM; ++i) {
    mean_values[i] /= (2.0 * number_lines);
  }
  if (number_lines == 3) {
    centroid = mean_values;
    centroid_calculated = true;
    return;
  }
  std::array<double, DIM> sum_up = {0};
  double sum_down = 0.0;
  for (int l = 0; l < number_lines; ++l) {
    const auto& line_start = lines[l]->get_start_point();
    const auto& line_end = lines[l]->get_end_point();
    for (int j = 0; j < DIM; ++j) {
      a[j] = line_start[j] - mean_values[j];
      b[j] = line_end[j] - mean_values[j];
      triangle_centroid[j] =
          (line_start[j] + line_end[j] + mean_values[j]) / 3.0;
    }
    const double A_l =
        0.5 * std::sqrt(std::pow(a[x2] * b[x3] - a[x3] * b[x2], 2.0) +
                        std::pow(a[x1] * b[x3] - a[x3] * b[x1], 2.0) +
                        std::pow(a[x2] * b[x1] - a[x1] * b[x2], 2.0));
    for (int i = 0; i < DIM; ++i) {
      sum_up[i] += A_l * triangle_centroid[i];
    }
    sum_down += A_l;
  }
  for (int i = 0; i < DIM; ++i) {
    centroid[i] = sum_up[i] / sum_down;
  }
  centroid_calculated = true;
}

void Polygon::calculate_normal() const {
  // Check if the centroid is calculated
  if (!centroid_calculated) {
    calculate_centroid();
  }
  // Find the normal vector for all the triangles formed from one edge of the
  // centroid
  for (int i = 0; i < number_lines; ++i) {
    for (int j = 0; j < DIM; ++j) {
      normals[i][j] = 0.0;
    }
  }
  std::array<double, DIM> v_out = {0};
  for (int i = 0; i < number_lines; ++i) {
    const auto& l1 = lines[i]->get_start_point();
    const auto& l2 = lines[i]->get_end_point();
    for (int j = 0; j < DIM; ++j) {
      a[j] = l1[j] - centroid[j];
      b[j] = l2[j] - centroid[j];
    }
    normals[i][x1] = 0.5 * (a[x2] * b[x3] - a[x3] * b[x2]);
    normals[i][x2] = -0.5 * (a[x1] * b[x3] - a[x3] * b[x1]);
    normals[i][x3] = 0.5 * (a[x1] * b[x2] - a[x2] * b[x1]);
    normals[i][const_i] = 0.0;
    const auto& o = lines[i]->get_outside_point();
    for (int j = 0; j < DIM; ++j) {
      v_out[j] = o[j] - centroid[j];
    }
    flip_normal_if_needed(normals[i], v_out);
  }
  for (int j = 0; j < DIM; ++j) {
    normal[j] = 0.0;
  }
  for (int i = 0; i < number_lines; ++i) {
    for (int j = 0; j < DIM; ++j) {
      normal[j] += normals[i][j];
    }
  }
  normal_calculated = true;
}

void Polygon::print(std::ofstream& file, std::array<double, DIM> position) {
  // Print the polygon to the file
  for (int i = 0; i < number_lines; i++) {
    const auto& p1 = lines[i]->get_start_point();
    const auto& p2 = lines[i]->get_end_point();
    file << position[x1] + p1[x1] << " " << position[x2] + p1[x2] << " "
         << position[x3] + p1[x3] << " " << position[x1] + p2[x1] << " "
         << position[x2] + p2[x2] << " " << position[x3] + p2[x3] << " "
         << position[x1] + centroid[x1] << " " << position[x2] + centroid[x2]
         << " " << position[x3] + centroid[x3] << std::endl;
  }
}

}  // namespace JetscapeCornelius
