//
//  JetClass.cpp
//  
//
//  Created by Abhijit Majumder on 10/6/16.
//
//

#include <iostream>
#include <complex>
#include <fstream>
#include <math.h>
#include <assert.h>
#include "JetClass.hpp"
#include "constants.h"
#include "JetScapeLogger.h"

/*
Jet::Jet(FourVector p)
{
    set_jet_p(p);
}
*/

VertexBase::~VertexBase()
{
  VERBOSESHOWER(9);
}

Parton::~Parton()
{
  VERBOSESHOWER(9);
}

//Make smarter later ...

Parton::Parton (int label, int id, int stat, double pt, double eta, double phi, double e)
{
  set_label(label);
    
  set_id(id);
  
  set_mass(-1.0);
  switch (id) {
  case 1:  //down quark
  case -1: // anti-down quark
    set_mass(0.01);
    break;
    
  case 2:  // up quark
  case -2:  // anti-up quark
    set_mass(0.005);
    break;
    
  case 3:   // strange quark
  case -3:  // anti-strange quark
    set_mass(0.15);
    break;
    
  case 21: // gluon
    set_mass(0.0);
    break;
    
  default:
    {
      std::cout << " error in id = " << id << std::endl;
      assert(mass()>=0.0);
      break;
    }
  }
  
    double p[4];
    p[0] = e;
    p[1] = pt*cos(phi);
    p[2] = pt*sin(phi);
    p[3] = pt*sinh(eta);
    
    
    set_p(p);

    set_stat(stat);
    
    double x[4];
    x[0]=0;
    x[1]=0;
    x[2]=0;
    x[3]=0;

    set_x(0); // if no x specified in constructor, particle starts at origin
  
  reset_PtYPhiM(pt,eta,phi,mass()); //check
}

Parton::Parton (int label, int id, int stat, double pt, double eta, double phi, double e, double x[4]) 
{
   set_label(label);
    
   set_id(id);
   
   set_mass(-1.0);
   switch (id)
    {

        case 1:  //down quark
        case -1: // anti-down quark
        set_mass(0.01);
        break;
     

        case 2:  // up quark
        case -2:  // anti-up quark
        set_mass(0.005);
        break;
    
        case 3:   // strange quark
        case -3:  // anti-strange quark
        set_mass(0.15);
        break;
     

        case 21: // gluon
        set_mass(0.0);
        break;
     
        default:
        {
            std::cout << " error in id = " << id << std::endl;
            assert(mass()>=0.0);
            break;
        }
    }
   
    double p[4];
    p[0] = e;
    p[1] = pt*cos(phi);
    p[2] = pt*sin(phi);
    p[3] = pt*sinh(eta);
    
    
    set_p(p);
    set_stat(stat);
    set_x(x);

   reset_PtYPhiM(pt,eta,phi,mass()); //check
}

Parton::Parton (int label, int id, int stat, double p[4], double x[4])  : PseudoJet(p[0],p[1],p[2],p[3])
{
    
    set_label(label);
    
    set_id(id);
    
    set_mass(-1.0);
    switch (id) {
        case 1:  //down quark
        case -1: // anti-down quark
            set_mass(0.01);
            break;
            
        case 2:  // up quark
        case -2:  // anti-up quark
            set_mass(0.005);
            break;
        
        case 3:   // strange quark
        case -3:  // anti-strange quark
            set_mass(0.15);
            break;
            
        case 21: // gluon
            set_mass(0.0);
            break;
            
        default:
        {
            std::cout << " error in id = " << id << std::endl;
            assert(mass()>=0.0);
            break;
        }
    }
    set_p(p);
    set_x(x);
    set_stat(stat);
}

