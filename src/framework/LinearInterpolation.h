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

#ifndef LINEARINTERPOLATION_H
#define LINEARINTERPOLATION_H

#include <cmath>

#include "RealType.h"

namespace Jetscape{

/// any type with + and scale * overloaded can use this function
template <class type>
type LinearInt(real x0, real x1, type y0, type y1, real x) {
    type temp = ((x - x0) * y1 + (x1 - x) * y0) / (x1 - x0);
    return temp;
}


// inspired by numerical recipes
// x0,x1: grid points in x-direction
// y0,y1: grid points in y-direction
// f0-f3: function value starting at x0,y0, continue counterclockwise
// put differently: f0=f(x0,y0)
// f1=f(x1,y0)
// f2=f(x1,y1)
// f3=f(x0,y1)
template <class type>
type BilinearInt(real x0, real x1, real y0, real y1,
                  type f0, type f1, type f2, type f3,
                  real x, real y) {
    type temp;
    real t = (x - x0)/(x1 - x0);
    real u = (y - y0)/(y1 - y0);
    if ((std::isfinite(u) == 1) && (std::isfinite(t) == 1)) {
        temp = (1 - t)*(1 - u)*f0+t*(1 - u)*f1 + t*u*f2 + (1 - t)*u*f3;
    } else {
        if (std::isfinite(u) == 0) temp = LinearInt(x0, x1, f0, f2, x);
        if (std::isfinite(t) == 0) temp = LinearInt(y0, y1, f0, f2, y);
    }
    return temp;
}


// 3D linear interpolation
template <class type>
type TrilinearInt(real x0, real x1, real y0, real y1, real z0, real z1,
                   type f000, type f001, type f010, type f011,
                   type f100, type f101, type f110, type f111,
                   real x, real y, real z) {
    type temp;
    real t = (x - x0)/(x1 - x0);
    real u = (y - y0)/(y1 - y0);
    real v = (z - z0)/(z1 - z0);

    if ((std::isfinite(u)==1)&&(std::isfinite(t)==1)&&(std::isfinite(v)==1)) {
        temp = (1-t)*(1-u)*(1-v)*f000;
        temp = temp + (1-t)*(1-u)*v*f001;
        temp = temp + (1-t)*u*(1-v)*f010;
        temp = temp + (1-t)*u*v*f011;
        temp = temp + t*(1-u)*(1-v)*f100;
        temp = temp + t*(1-u)*v*f101;
        temp = temp + t*u*(1-v)*f110;
        temp = temp + t*u*v*f111;
    } else {
        if (std::isfinite(t)==0)
            temp=BilinearInt(y0,y1,z0,z1,f000,f010,f011,f001,y,z);
        if (std::isfinite(v)==0)
            temp=BilinearInt(x0,x1,y0,y1,f000,f100,f110,f010,x,y);
        if (std::isfinite(u)==0)
            temp=BilinearInt(x0,x1,z0,z1,f000,f100,f101,f001,x,z);
    }
    return temp;
}

} //end namespace Jetscape

#endif  // LINEARINTERPOLATION_H
