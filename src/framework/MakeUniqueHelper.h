/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
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
 * Allow C++14 style make_unique for C++11
 * Taken from STL's solution,
 * https://isocpp.org/files/papers/N3656.txt
 ******************************************************************************/


#ifndef MAKEUNIQUEHELPER_H
#define MAKEUNIQUEHELPER_H

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

namespace Jetscape {
  template<class T> struct _Unique_if {
    typedef std::unique_ptr<T> _Single_object;
  };

  template<class T> struct _Unique_if<T[]> {
    typedef std::unique_ptr<T[]> _Unknown_bound;
  };

  template<class T, size_t N> struct _Unique_if<T[N]> {
    typedef void _Known_bound;
  };

  template<class T, class... Args>
    typename _Unique_if<T>::_Single_object
    make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

  template<class T>
    typename _Unique_if<T>::_Unknown_bound
    make_unique(size_t n) {
    typedef typename std::remove_extent<T>::type U;
    return std::unique_ptr<T>(new U[n]());
  }

  template<class T, class... Args>
    typename _Unique_if<T>::_Known_bound
    make_unique(Args&&...) = delete;

  // check whether a weak pointer is initialized or not
  template <typename T>
    bool weak_ptr_is_uninitialized(std::weak_ptr<T> const& weak) {
    using wt = std::weak_ptr<T>;
    return !weak.owner_before(wt{}) && !wt{}.owner_before(weak);
  }

}

#endif // MAKEUNIQUEHELPER_H
