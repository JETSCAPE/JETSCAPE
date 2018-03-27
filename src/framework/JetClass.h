//
//  JetClass.h
//  
//
//  Created by Abhijit Majumder on 10/6/16.
//
//

#ifndef JetClass_h
#define JetClass_h

#include <stdio.h>
#include <math.h>
#include "JetScapeParticles.hpp"
#include "JetScapeConstants.h"
#include "four_vector.hpp"
#include "fjcore.hh"

#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>

using std::ostream;

namespace Jetscape {


  // class Parton;
class Vertex;
class FourVector;


/**************************************************************************************************/

//  JET CLASS

/*************************************************************************************************/

//dummy for now figure out after graph structure ...

class Jet
{
  Jet() {};
  ~Jet() {};
};


/**************************************************************************************************/

//  VERTEX CLASS

/*************************************************************************************************/

class Vertex
{
    
public:

  Vertex() {x_in_.Set(0,0,0,0);}
  Vertex(double x, double y, double z, double t) {x_in_.Set(x,y,z,t);}
  Vertex(FourVector &x) {set_location(x);}  
  virtual ~Vertex();
  
  void set_location(FourVector &x)
    {
      x_in_ = x;
    }
  
  FourVector &x_in()
  {
    return(x_in_);
  }

  friend ostream &operator<<( ostream &output, 
         Vertex & vertex ) {

    output<<vertex.x_in().x()<<" "<<vertex.x_in().y()<<" "<<vertex.x_in().z()<<" "<<vertex.x_in().t();//<<endl;
    
    return output;
  }
protected:
  
  FourVector x_in_        ; //location of the vertex
  // parents and siblings from Graph structure later ...
  
};

};  /// end of namespace Jetscape

#endif /* JetClass_h */
