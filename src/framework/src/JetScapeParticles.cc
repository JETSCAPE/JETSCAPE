/** Provides JetScapeParticleBase and derived classes Parton, Hadron

    @group JetScape (modular/task) based framework
    @author Kolja Kauder
    @version Revision 0.1
    @date Aug 8, 2017
*/

#include <iostream>
#include <complex>
#include <fstream>
#include <cmath>
#include <assert.h>
#include "JetScapeParticles.hpp"
#include "constants.h"
#include "JetScapeLogger.h"

namespace Jetscape {

  JetScapeParticleBase::~JetScapeParticleBase(){
    VERBOSESHOWER(9);
  }

  JetScapeParticleBase::JetScapeParticleBase (const JetScapeParticleBase& srp)
    : PseudoJet (srp)
  {
    pid_ = srp.pid_;
    plabel_ = srp.plabel_;
    pstat_ = srp.pstat_;
    form_time_ = srp.form_time_;
    mass_ = srp.mass_;
    // t_ = srp.t_;
    jet_v_ = srp.jet_v_;
    x_in_ = srp.x_in_;
  }

  JetScapeParticleBase::JetScapeParticleBase (int label, int id, int stat, double pt, double eta, double phi, double e, double* x)
  {
    set_label(label);
    set_id(id);
    initialize_form_time();
    init_jet_v();
  
    set_restmass(-1.0);
    switch (id) {
    case 1:  //down quark
    case -1: // anti-down quark
      set_restmass(0.01);
      break;
    
    case 2:  // up quark
    case -2:  // anti-up quark
      set_restmass(0.005);
      break;
    
    case 3:   // strange quark
    case -3:  // anti-strange quark
      set_restmass(0.15);
      break;
    
    case 4:   // charm quark
    case -4:  // anti-charm quark
      set_restmass(1.29); // make more accurate later
      break;

    case 5:   // bottom quark
    case -5:  // anti-bottom quark
      set_restmass(4.2); // make more accurate later
      break;
      
    case 21: // gluon
      set_restmass(0.0);
      break;
      
    default:
      std::cout << " error in id = " << id << std::endl;
      assert(mass_>=0.0);
      break;
    }
    
    double p[4];
    p[0] = e;
    p[1] = pt*cos(phi);
    p[2] = pt*sin(phi);
    p[3] = pt*sinh(eta);
    
    set_p(p);
    set_stat(stat);

    if ( x ){
      set_x(x);
    } else {
      // if no x specified in constructor, particle starts at origin
      double x0[4];
      x0[0]=0;
      x0[1]=0;
      x0[2]=0;
      x0[3]=0;
      set_x(x0); 
    }
  
    // reset_PtYPhiM(pt,eta,phi,mass_); //check
  }

    
  JetScapeParticleBase::JetScapeParticleBase (int label, int id, int stat, double p[4], double x[4])  : PseudoJet(p[0],p[1],p[2],p[3])
  {
    
    set_label(label);
    set_id(id);    
    initialize_form_time();    
    init_jet_v();
    
    set_restmass(-1.0);
    switch (id) {
    case 1:  //down quark
    case -1: // anti-down quark
      set_restmass(0.01);
      break;
            
    case 2:  // up quark
    case -2:  // anti-up quark
      set_restmass(0.005);
      break;
        
    case 3:   // strange quark
    case -3:  // anti-strange quark
      set_restmass(0.15);
      break;

    case 4:   // charm quark
    case -4:  // anti-charm quark
      set_restmass(1.29); // make more accurate later
      break;
            
    case 5:   // bottom quark
    case -5:  // anti-bottom quark
      set_restmass(4.2); // make more accurate later
      break;
            

    case 21: // gluon
      set_restmass(0.0);
      break;
            
    default:
      {
	std::cout << " error in id = " << id << std::endl;
	assert(mass_>=0.0);
	break;
      }
    }
    set_p(p);
    set_x(x);
    set_stat(stat);
  }


  void JetScapeParticleBase::clear()
  {
    plabel_ = 0;
    pid_ = 0;
    pstat_ = 0;        
  }
    
  void JetScapeParticleBase::set_label(int label)
  {
    plabel_ = label;
  }
  
  void JetScapeParticleBase::set_id(int id)
  {
    pid_ = id;
  }
  
  void JetScapeParticleBase::set_stat(int stat)
  {
    pstat_ = stat;
  }

  void JetScapeParticleBase::set_restmass(double mass_input)
  {
    mass_ = mass_input;
  }
    
  void JetScapeParticleBase::set_p(double p[4])
  {
    reset_momentum(p[0],p[1],p[2],p[3]);
    // p_in_ = FourVector(p);
    // //FourVector p_in_(p); // error: creates new vector and hence not accessible via class p_in_
    // //reset_momentum(p[0],p[1],p[2],p[3]);
  }

  // not needed in graph structure
  
  void JetScapeParticleBase::set_x(double x[4])
  {
    //FourVector
    x_in_.Set(x);
  }
 
  void JetScapeParticleBase::set_mean_form_time ()
  {
    //mean_form_time_ = (this->e()+this->pl())*std::sqrt(2.0)/t_;
    // mean_form_time_ = 2.0*this->e()/t_;
    mean_form_time_ = 2.0*e()/t();
  }
    

  void JetScapeParticleBase::set_form_time(double form_time)
  {
    form_time_ = form_time;
  }
    
  void JetScapeParticleBase::initialize_form_time()
  {
    form_time_ = -0.1;
  }
    
  void JetScapeParticleBase::set_t(double t)
  {
    //  Reset the momentum due to virtuality              
    double newP[4];              
    newP[0] = e();
    double newPl = std::sqrt( e()*e() - t ) ;              
    double velocityMod = std::sqrt(std::pow(jet_v_.comp(1),2) + std::pow(jet_v_.comp(2),2) + std::pow(jet_v_.comp(3),2));

    newPl = newPl/velocityMod;
    for(int j=1;j<=3;j++) {
      newP[j] = newPl*jet_v_.comp(j);
    }
    
    set_p(newP);

	      
    // t_ = t;
  } 

    
  void JetScapeParticleBase::init_jet_v()
  {
    jet_v_ = FourVector();
  }
    
  void JetScapeParticleBase::set_jet_v(double v[4])
  {
    jet_v_ =FourVector(v);
  }
    
  //  end setters
    
  //  start getters
    
  double JetScapeParticleBase::form_time()
  {
    return(form_time_);
  }
    
  const int JetScapeParticleBase::pid()
  {
    return(pid_);
  }
    
  const int JetScapeParticleBase::pstat()
  {
    return(pstat_);
  }
    
  const int JetScapeParticleBase::plabel()
  {
    return(plabel_);
  }

  // const double JetScapeParticleBase::e()
  // {
  //   return(p_in_.t());
  // }
    
  // const double JetScapeParticleBase::pt()
  // {
  //   return(sqrt(p_in_.x()*p_in_.x() + p_in_.y()*p_in_.y())) ;    
  // }
    
    
  const double JetScapeParticleBase::time()
  {
    return(x_in_.t());
  }
    
  /*
    int pparent_label()
    {
    return(pparent_label_);
    }
  */
  // FourVector &JetScapeParticleBase::p_in()
  // {
  //   return(p_in_);
  // }
  
  
  FourVector &JetScapeParticleBase::x_in()
  {
    return(x_in_);
  }
    
  FourVector &JetScapeParticleBase::jet_v()
  {
    return(jet_v_);
  }
  
  const double JetScapeParticleBase::restmass()
  {
    return(mass_);
  }

  const double JetScapeParticleBase::mean_form_time()
  {
    return(mean_form_time_);
  }

  // just operator of PseudoJet ...
  
  const double JetScapeParticleBase::p(int i) {
    // return (p_in_.comp(i));
    return operator()(i);
  }
  
  
  double JetScapeParticleBase::pl() {
    if (jet_v_.comp(0)<0.99) {
      // this should never happen
      cerr << "jet_v_ = " << jet_v_.comp(0) << "  "  << jet_v_.comp(1) << "  "  << jet_v_.comp(2) << "  "  << jet_v_.comp(3) << endl;
      throw std::runtime_error("JetScapeParticleBase::pl() : jet_v should never be space-like.");
      // return(std::sqrt( p_in_.x()*p_in_.x() + p_in_.y()*p_in_.y() + p_in_.z()*p_in_.z() ) );
      return(-1);
    }
    else
      {
	// cout << " returing pl using jet_v " << ( p_in_.x()*jet_v_.x()+ p_in_.y()*jet_v_.y() + p_in_.z()*jet_v_.z() )/std::sqrt( pow(jet_v_.x(),2) + pow(jet_v_.y(),2) + pow(jet_v_.z(),2) ) << endl ;
	// return( (p_in_.x()*jet_v_.x()+ p_in_.y()*jet_v_.y() + p_in_.z()*jet_v_.z() )/std::sqrt( pow(jet_v_.x(),2) + pow(jet_v_.y(),2) + pow(jet_v_.z(),2) ) );
	// projection onto (unit) jet velocity
	return ( px()*jet_v_.x()+ py()*jet_v_.y() + pz()*jet_v_.z() )/std::sqrt( pow(jet_v_.x(),2) + pow(jet_v_.y(),2) + pow(jet_v_.z(),2) );
      }
  }
    
  const double JetScapeParticleBase::nu()
  {
    return( ( this->e()+std::abs(this->pl()) )/std::sqrt(2) );
  }
    
  const double JetScapeParticleBase::t()
  {
    /// \Todo: Fix 
    return ( PseudoJet::m2() ) ;
    // return (t_) ;
  }        
    
  /*   double generate_t(double, double);
       {
    
       return (1);
       }
  */
 
  JetScapeParticleBase& JetScapeParticleBase::operator=(JetScapeParticleBase &c)
  {
    fjcore::PseudoJet::operator=(c);
    //FourVector x_in_;
      
    pid_ = c.pid() ;
    pstat_ = c.pstat() ;
    plabel_ = c.plabel() ;
    //pparent_label_ = c.pparent_label();
    // p_in_ = c.p_in() ;
    //JetScapeParticleBase::set_x(c.x_in());

    x_in_ = c.x_in() ;
    form_time_ = c.form_time_;
    mass_ = c.mass_;
    // t_ = c.t();
      
    return *this;
  }
  
  JetScapeParticleBase& JetScapeParticleBase::operator=(const JetScapeParticleBase &c)
  {
    fjcore::PseudoJet::operator=(c);
    //FourVector x_in_;
      
    pid_ = c.pid_;
    pstat_ = c.pstat_ ;
    plabel_ = c.plabel_;
    //pparent_label_ = c.pparent_label();
    // p_in_ = c.p_in_;
    //  x_in_ = c.x_in() ;
        
    //JetScapeParticleBase::set_x(c.x_in());

    x_in_ = c.x_in_;
      
    mass_ = c.mass_;
    // t_ = c.t_;
    form_time_ = c.form_time_ ;
        
    return *this;
  }
  /*
    friend ostream Print()
    {
    ostream output;
    return output;
    }
  */
  
  ostream &operator<<( ostream &output, JetScapeParticleBase & p ) {
    output<<p.plabel()<<" "<<p.pid()<<" "<<p.pstat()<<" ";
    output<<p.pt()<<" "<<p.rap()<<" "<<p.phi()<<" "<<p.e()<<" ";
    output<<p.x_in().x()<<" "<<p.x_in().y()<<" "<<p.x_in().z()<<" "<<p.x_in().t();//<<endl;
    
    return output;            
  }

} /// end of namespace Jetscape
