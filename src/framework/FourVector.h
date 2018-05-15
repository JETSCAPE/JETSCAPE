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

#ifndef FOURVECTOR_H
#define FOURVECTOR_H

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <climits>

using std::cout;
using std::cerr;
using std::endl;

namespace Jetscape {

  class FourVector
  {
    // the class of four vectors



  public:
    
    
    FourVector() //default constructor
    {
      tv=xv=yv=zv=0.0;
    };
    
    FourVector(const FourVector& srv): xv(srv.xv) , yv(srv.yv), zv(srv.zv), tv(srv.tv)
    {};  // copy constructor
    
    FourVector(double a[4])  // constructor with array input
    {
      tv=a[0];
      xv=a[1];
      yv=a[2];
      zv=a[3];
    };

    
    FourVector(double x_in, double y_in, double z_in, double t_in) // constructor with component input
    {
      tv=t_in;
      xv=x_in;
      yv=y_in;
      zv=z_in;
    };
    
    void clear()
    {
      tv=xv=yv=zv=0.0;
    }
    
    
    // constructors do all sets

    void Set(double x_in, double y_in, double z_in, double t_in)
    {
      tv=t_in;
      xv=x_in;
      yv=y_in;
      zv=z_in;
    }
  
    void Set(double a[4]) 
    {
      tv=a[0];
      xv=a[1];
      yv=a[2];
      zv=a[3];
    };
    
    // all gets are done with name calls e.g., vec.x()
    double x() const
    {
      return(xv);
    };
    
    double y() const
    {
      return(yv);
    };
    
    double z() const
    {
      return(zv);
    };
    
    double t() const
    {
      return(tv);
    };
    
    const double comp(int i) const
    {
      switch (i) {
      case 0:
	return(tv);
	break;
      case 1:
	return(xv);
	break;
      case 2:
	return(yv);
	break;
      case 3:
	return(zv);
	break;
      default:
	cout << " component index beyond 0-3! Returning garbage ..." << endl ;
	return(a_very_large_number);
	break;
      }
    }
  
    double plus()
    {
      return ( (zv+tv)/sqrt(2.0) );
    };
    
    double minus()
    {
      return ( (tv-zv)/sqrt(2.0) );
    };
    
    double rapidity()
    {
      if (this->minus()>0.0) return ( std::log(this->plus()/this->minus() )/2.0  );
      cout << endl << "ERROR: z component exceeds t component, cannot calculate rapidity" << endl;
      return (0);
    };
    
    double phi() {
      if ( fabs(x())<rounding_error && fabs(y())<rounding_error ) {
	return 0;
      }
      double phi=atan2( y(), x() );
      while (phi<0) phi+=2.0*pi;
      return phi;
    };

    double operator*(FourVector &c)
    {
      return(tv*c.t() - xv*c.x() - yv*c.y() - zv*c.z() );
    };
    
    
    FourVector &operator+=(FourVector &c)
    {
      tv+=c.t();
      xv+=c.x();
      yv+=c.y();
      zv+=c.z();
        
      return(*this);
    };

    FourVector &operator-=(FourVector &c)
    {
        tv-=c.t();
        xv-=c.x();
        yv-=c.y();
        zv-=c.z();
        
        return(*this);
    };

    FourVector &operator=(FourVector &c)
    {
      tv = c.t();
      xv = c.x();
      yv = c.y();
      zv = c.z();
      return (*this);
        
    };

    FourVector &operator=(const FourVector &c)
    {
      tv = c.tv;
      xv = c.xv;
      yv = c.yv;
      zv = c.zv;
      return (*this);
        
    };
    
    void rotate_around_z(double theta)
    {
      double new_xv, new_yv;
        
      new_xv = xv*cos(theta) - yv*sin(theta);
        
      new_yv = yv*cos(theta) + xv*sin(theta);
        
      xv = new_xv;
      yv = new_yv;
    };

  private:
    // the v is for vector, we call the private variables, xv, tv etc., so that get function
    // calls will be called x, t etc.
    double xv,yv,zv,tv;
   

  };


}; /// end of namespace Jetscape

#endif // FOURVECTOR_H
