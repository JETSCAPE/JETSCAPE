/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

/*******************************************************************************
 * @file MakeUniqueHelper.h
 * @brief Helper utilities for C++11 to emulate C++14's std::make_unique.
 *
 * This header defines a `make_unique` implementation that provides equivalent functionality
 * to C++14's `std::make_unique`, allowing safer dynamic memory allocation with `std::unique_ptr`
 * in C++11. It includes support for single objects and dynamically sized arrays, but prevents
 * use with fixed-size arrays, as per the C++14 standard proposal (N3656).
 *
 * It also provides a utility function to determine if a `std::weak_ptr` is uninitialized.
 *
 * Adapted from: https://isocpp.org/files/papers/N3656.txt
 *
 * @copyright
 * Copyright (c) The JETSCAPE Collaboration, 2018
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 ******************************************************************************/

#ifndef MAKEUNIQUEHELPER_H
#define MAKEUNIQUEHELPER_H

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

namespace Jetscape {

/**
 * @brief Internal trait to determine the correct return type for make_unique.
 *
 * Specializes for single objects, dynamically sized arrays, and fixed-size arrays.
 */
template <class T>
struct _Unique_if {
  typedef std::unique_ptr<T> _Single_object; ///< For single object allocation
};

template <class T>
struct _Unique_if<T[]> {
  typedef std::unique_ptr<T[]> _Unknown_bound; ///< For dynamically sized arrays
};

template <class T, size_t N>
struct _Unique_if<T[N]> {
  typedef void _Known_bound; ///< Disallowed: fixed-size arrays
};

/**
 * @brief Create a unique_ptr to a single object.
 *
 * @tparam T Object type
 * @tparam Args Constructor argument types
 * @param args Arguments to forward to T's constructor
 * @return std::unique_ptr<T>
 */
template <class T, class... Args>
typename _Unique_if<T>::_Single_object make_unique(Args &&...args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

/**
 * @brief Create a unique_ptr to a dynamically allocated array.
 *
 * @tparam T Array element type
 * @param n Number of elements in the array
 * @return std::unique_ptr<T[]> to a zero-initialized array
 */
template <class T>
typename _Unique_if<T>::_Unknown_bound make_unique(size_t n) {
  typedef typename std::remove_extent<T>::type U;
  return std::unique_ptr<T>(new U[n]());
}

/**
 * @brief Deleted function to disallow make_unique for fixed-size arrays.
 */
template <class T, class... Args>
typename _Unique_if<T>::_Known_bound make_unique(Args &&...) = delete;

/**
 * @brief Utility to check if a std::weak_ptr is uninitialized.
 *
 * @tparam T The type pointed to
 * @param weak The weak pointer to check
 * @return true if the weak pointer has not been initialized, false otherwise
 */
template <typename T>
bool weak_ptr_is_uninitialized(std::weak_ptr<T> const &weak) {
  using wt = std::weak_ptr<T>;
  return !weak.owner_before(wt{}) && !wt{}.owner_before(weak);
}

}  // namespace Jetscape

#endif  // MAKEUNIQUEHELPER_H
