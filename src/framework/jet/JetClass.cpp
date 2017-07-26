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
#include <cmath>
#include <assert.h>
#include "JetClass.hpp"
#include "constants.h"
#include "JetScapeLogger.h"

namespace Jetscape {



/*
Jet::Jet(FourVector p)
{
    set_jet_p(p);
}
*/

Vertex::~Vertex()
{
  VERBOSESHOWER(9);
}

Parton::~Parton()
{
  VERBOSESHOWER(9);
}

//Make smarter later ...


Parton::Parton (const Parton& srp)
{
    pid_ = srp.pid_;
    pstat_ = srp.pstat_;
    form_time_ = srp.form_time_;
    mass_ = srp.mass_;
    t_ = srp.t_;
    t_max_ = srp.t_max_;
    jet_v_ = srp.jet_v_;
    p_in_ = srp.p_in_;
    x_in_ = srp.x_in_;
//    cout << " form time in cpC = " << form_time_ << endl ;
};

Parton::Parton (int label, int id, int stat, double pt, double eta, double phi, double e)
{
  set_label(label);
    
  set_id(id);

  initialize_form_time();
    
    init_jet_v();
  
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
    
  case 4:   // charm quark
  case -4:  // anti-charm quark
    set_mass(1.29); // make more accurate later
    break;

  case 5:   // bottom quark
  case -5:  // anti-bottom quark
    set_mass(4.2); // make more accurate later
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

    set_x(x); // if no x specified in constructor, particle starts at origin
    
  
  reset_PtYPhiM(pt,eta,phi,mass()); //check
}

Parton::Parton (int label, int id, int stat, double pt, double eta, double phi, double e, double x[4]) 
{
   set_label(label);
    
   set_id(id);
    
   initialize_form_time();

    init_jet_v();

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
     
        case 4:   // charm quark
        case -4:  // anti-charm quark
            set_mass(1.29); // make more accurate later
            break;
            
        case 5:   // bottom quark
        case -5:  // anti-bottom quark
            set_mass(4.2); // make more accurate later
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
    
    initialize_form_time();
    
    init_jet_v();

    
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

        case 4:   // charm quark
        case -4:  // anti-charm quark
            set_mass(1.29); // make more accurate later
            break;
            
        case 5:   // bottom quark
        case -5:  // anti-bottom quark
            set_mass(4.2); // make more accurate later
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


/*double Parton::generate_t(double t_min, double t_max)
{
    
   // cout << " inside generate parton momentum " << this->p_in_.t() << "  " << this->p_in_.x() << "  " << this->p_in_.y() << "  " << this->p_in_.z() << endl ;
   // cout << std::pow(this->p_in_.t(),2) << "  "  << std::pow(this->p_in_.x(),2) << "  " << std::pow(this->p_in_.y(),2) << "  " << std::pow(this->p_in_.z(),2) << endl ;
    
   double tPrior = std::pow(this->p_in_.t(),2) - std::pow(this->p_in_.x(),2) - std::pow(this->p_in_.y(),2) - std::pow(this->p_in_.z(),2) ;
    
    tPrior = 1.0;
    
    double loc = (this->x_in_*this->p_in_)/this->p(0);
    
   // cout << " tPrior = " << tPrior << endl;
    if (t_<-99.0)
    {
        
        double nu = (this->e() + this->pl())/std::sqrt(2);
        t_ = generate_vac_t(this->pid_,nu,t_min, t_max,loc)
        
    }
    
    return (tPrior);
    
}

*/


} /// end of namespace Jetscape
