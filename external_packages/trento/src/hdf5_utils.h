// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// TRENTO3D: Three-dimensional extension of TRENTO by Weiyao Ke
// MIT License

#ifndef HDF5_UTILS_H
#define HDF5_UTILS_H

#include <string>

#ifdef TRENTO_HDF5
#include <H5Cpp.h>
// This macro was introduced in v1.8.7.  Define it manually for older versions.
#ifndef H5_VERSION_GE
#define H5_VERSION_GE(Maj,Min,Rel) \
       (((H5_VERS_MAJOR==Maj) && (H5_VERS_MINOR==Min) && (H5_VERS_RELEASE>=Rel)) || \
        ((H5_VERS_MAJOR==Maj) && (H5_VERS_MINOR>Min)) || \
        (H5_VERS_MAJOR>Maj))
#endif
#endif

#include "fwd_decl.h"

namespace trento {

namespace hdf5 {

// Determine if a filename is an HDF5 file based on the extension.
bool filename_is_hdf5(const fs::path& path);
bool filename_is_hdf5(const std::string& path);

#ifdef TRENTO_HDF5

// Open an HDF5 file object.
// Throw std::invalid_argument if the file does not exist or is not valid HDF5.
H5::H5File try_open_file(
    const std::string& path, unsigned int flags = H5F_ACC_RDONLY);

// Map C types to corresponding HDF5 datatypes.
// See section "predefined datatypes" in the HDF5 docs.
using H5::PredType;
template <typename T> inline const PredType& type();
template <> inline const PredType& type<int>()           { return PredType::NATIVE_INT; }
template <> inline const PredType& type<unsigned long>() { return PredType::NATIVE_UINT; }
template <> inline const PredType& type<long int>()      { return PredType::NATIVE_LONG; }
template <> inline const PredType& type<long long int>() { return PredType::NATIVE_LLONG; }
template <> inline const PredType& type<float>()         { return PredType::NATIVE_FLOAT; }
template <> inline const PredType& type<double>()        { return PredType::NATIVE_DOUBLE; }
template <> inline const PredType& type<long double>()   { return PredType::NATIVE_LDOUBLE; }

// Construct an HDF5 "simple" dataspace from a generic container.
// It must contain values of type "hsize_t" to match the H5::DataSpace ctor.
template <typename Container>
inline typename std::enable_if<
  std::is_same<typename Container::value_type, hsize_t>::value,
  H5::DataSpace
>::type
make_dataspace(const Container& shape) {
  return H5::DataSpace{shape.size(), shape.data()};
}

#endif  // TRENTO_HDF5

}  // namespace hdf5

}  // namespace trento

#endif  // HDF5_UTILS_H
