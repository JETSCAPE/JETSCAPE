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

#ifndef JetClass_h
#define JetClass_h

#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <vector>

#include "FourVector.h"
#include "JetScapeConstants.h"
#include "JetScapeParticles.h"
#include "fjcore.hh"

using std::ostream;

namespace Jetscape {

// class Parton;
class Vertex;
class FourVector;

/**************************************************************************************************/

//  JET CLASS

/*************************************************************************************************/

// dummy for now figure out after graph structure ...

class Jet {
  Jet(){};
  ~Jet(){};
};

/**************************************************************************************************/

//  VERTEX CLASS

/*************************************************************************************************/

class Vertex {
 public:
  Vertex() { x_in_.Set(0, 0, 0, 0); }
  Vertex(double x, double y, double z, double t) { x_in_.Set(x, y, z, t); }
  Vertex(FourVector &x) { set_location(x); }
  virtual ~Vertex();

  void set_location(FourVector &x) { x_in_ = x; }

  FourVector &x_in() { return (x_in_); }

  friend ostream &operator<<(ostream &output, Vertex &vertex) {
    output << vertex.x_in().x() << " " << vertex.x_in().y() << " "
           << vertex.x_in().z() << " " << vertex.x_in().t();
    return output;
  }

 protected:
  FourVector x_in_;  // location of the vertex
  // parents and siblings from Graph structure later ...
};

};  // namespace Jetscape

#endif /* JetClass_h */
