//
//  plotter.h
//  
//
//  Created by Abhijit Majumder on 7/28/16.
//
//

#ifndef _plotter_h
#define _plotter_h


#include <complex>
#include <math.h>
#include <stdlib.h>

namespace Jetscape {

class sudakov
{
public:
    double value (double cg0, double cg);
    double value_ex (double cg0, double cg);
};


} // end namespace Jetscape
#endif
