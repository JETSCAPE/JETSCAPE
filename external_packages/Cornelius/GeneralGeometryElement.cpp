#include "GeneralGeometryElement.h"

namespace JetscapeCornelius {

GeneralGeometryElement::GeneralGeometryElement()
    : normal_calculated(false), centroid_calculated(false) {}

GeneralGeometryElement::~GeneralGeometryElement() = default;

void GeneralGeometryElement::calculate_normal() const {
  // Provide implementation in derived classes
}

void GeneralGeometryElement::calculate_centroid() const {
  // Provide implementation in derived classes
}

}  // namespace JetscapeCornelius
