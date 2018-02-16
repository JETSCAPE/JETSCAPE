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
#include "JetScapeLogger.h"
#include "JetScapeParticles.hpp"
#include "constants.h"

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
    mass_ = srp.mass_;
    // t_ = srp.t_;
    jet_v_ = srp.jet_v_;
    x_in_ = srp.x_in_;
  }

  JetScapeParticleBase::JetScapeParticleBase (int label, int id, int stat, double pt, double eta, double phi, double e, double* x)
  {
    set_label(label);
    set_id(id);
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
      std::cerr << " error in id = " << id << std::endl;
      assert(mass_>=0.0);
      break;
    }
    
    // double p[4];
    // p[0] = e;
    // p[1] = pt*cos(phi);
    // p[2] = pt*sin(phi);
    // p[3] = pt*sinh(eta);
    // reset_momentum( FourVector ( p ) ); // kk: also works
    reset_momentum( pt*cos(phi),pt*sin(phi), pt*sinh(eta), e );
    
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

    
  JetScapeParticleBase::JetScapeParticleBase (int label, int id, int stat, const FourVector& p, const FourVector& x)
  {
    
    set_label(label);
    set_id(id);    
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

    reset_momentum(p);
    x_in_=x;
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
    
  // not needed in graph structure
  
  void JetScapeParticleBase::set_x(double x[4])
  {
    //FourVector
    x_in_.Set(x);
  }
    
  void JetScapeParticleBase::init_jet_v()
  {
    jet_v_ = FourVector();
  }
    
  void JetScapeParticleBase::set_jet_v(double v[4])
  {
    jet_v_ =FourVector(v);
  }
    
    void JetScapeParticleBase::set_jet_v(FourVector j)
    {
        jet_v_ = j;
    }
    
  //  end setters
    
  //  start getters        
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

  // FourVector JetScapeParticleBase::get_p() const{
  //   return FourVector ( px(), py(), pz(), e() );
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

  // just operator of PseudoJet ...
  
  const double JetScapeParticleBase::p(int i) {
    // return (p_in_.comp(i));    
    // return operator()(i);
    // cerr << " DON'T USE ME VERY OFTEN!!" << endl;
    switch ( i ){
    case 0 :      return e();
    case 1 :      return px();
    case 2 :      return py();
    case 3 :      return pz();
    default :     throw std::runtime_error("JetScapeParticleBase::p(int i) : i is out of bounds.");
    }
  }
  
  
  double JetScapeParticleBase::pl() {
    // Have to catch the ones initialized to 0
    if (jet_v_.comp(0)<1e-6) {
      return(std::sqrt( px()*px() + py()*py() + pz()*pz() ) );
    }
    
    if (jet_v_.comp(0)<0.99) {
      // this should never happen
      cerr << "jet_v_ = " << jet_v_.comp(0) << "  "  << jet_v_.comp(1) << "  "  << jet_v_.comp(2) << "  "  << jet_v_.comp(3) << endl;
      throw std::runtime_error("JetScapeParticleBase::pl() : jet_v should never be space-like.");
      // return(std::sqrt( p_in_.x()*p_in_.x() + p_in_.y()*p_in_.y() + p_in_.z()*p_in_.z() ) );
      return(-1);
    } else {
      // projection onto (unit) jet velocity
      return ( px()*jet_v_.x()+ py()*jet_v_.y() + pz()*jet_v_.z() )/std::sqrt( pow(jet_v_.x(),2) + pow(jet_v_.y(),2) + pow(jet_v_.z(),2) );
    }
  }
  
  const double JetScapeParticleBase::nu()
  {
    return( ( this->e()+std::abs(this->pl()) )/std::sqrt(2) );
  }
        
  JetScapeParticleBase& JetScapeParticleBase::operator=(JetScapeParticleBase &c)
  {
    fjcore::PseudoJet::operator=(c);
      
    pid_ = c.pid() ;
    pstat_ = c.pstat() ;
    plabel_ = c.plabel() ;

    x_in_ = c.x_in() ;
    mass_ = c.mass_;
      
    return *this;
  }
  
  JetScapeParticleBase& JetScapeParticleBase::operator=(const JetScapeParticleBase &c)
  {
    fjcore::PseudoJet::operator=(c);
      
    pid_ = c.pid_;
    pstat_ = c.pstat_ ;
    plabel_ = c.plabel_;

    x_in_ = c.x_in_;
      
    mass_ = c.mass_;        
    return *this;
  }
  
  ostream &operator<<( ostream &output, JetScapeParticleBase & p ) {
    output<<p.plabel()<<" "<<p.pid()<<" "<<p.pstat()<<" ";
    // output<<p.pt()<<" "<< (fabs (p.rap())>1e-15?p.rap():0)<<" "<< p.phi() <<" "<<p.e()<<" ";
    output<<p.pt()<<" "<< (fabs (p.eta())>1e-15?p.eta():0)<<" "<< p.phi() <<" "<<p.e()<<" ";
    output<<p.x_in().x()<<" "<<p.x_in().y()<<" "<<p.x_in().z()<<" "<<p.x_in().t();//<<endl;
    
    return output;            
  }

  // ---------------
  // Parton specific
  // ---------------
  
  Parton::Parton (const Parton& srp) :
    JetScapeParticleBase::JetScapeParticleBase (srp)
  {
    form_time_ = srp.form_time_;
      Color_ = srp.Color_;
      antiColor_ = srp.antiColor_;
      MaxColor_ = srp.MaxColor_;
      MinColor_ = srp.MinColor_;
      MinAntiColor_ = srp.MinAntiColor_;
      
    set_edgeid ( srp.edgeid() );
    pShower_ = srp.shower();
    // set_edgeid( -1 ); // by default do NOT copy the shower or my position in it
    // pShower_ = nullptr;
  }

  Parton::Parton (int label, int id, int stat, const FourVector& p, const FourVector& x)  :
    JetScapeParticleBase::JetScapeParticleBase ( label,  id,  stat,  p, x)
  {
    initialize_form_time();
      set_color(0);
      set_anti_color(0);
      set_min_color(0);
      set_min_anti_color(0);
      set_max_color(0);
    set_edgeid( -1 );
    pShower_ = nullptr;
    //   cout << "========================== std Ctor called, returning : " << endl << *this << endl;
  }
  

  Parton::Parton (int label, int id, int stat, double pt, double eta, double phi, double e, double* x)  :
    JetScapeParticleBase::JetScapeParticleBase ( label,  id,  stat,  pt, eta, phi, e, x){
    initialize_form_time();
        set_color(0);
        set_anti_color(0);
        set_min_color(0);
        set_min_anti_color(0);
        set_max_color(0);
    set_edgeid( -1 );
    pShower_ = nullptr;
    // cout << "========================== phieta Ctor called, returning : " << endl << *this << endl;
  }
  
  Parton& Parton::operator=( Parton &c)
  {
    JetScapeParticleBase::operator=(c);
    form_time_ = c.form_time_;
      Color_ = c.Color_;
      antiColor_ = c.antiColor_;
    set_edgeid ( c.edgeid() );
    pShower_ = c.shower();
    // set_edgeid( -1 ); // by default do NOT copy the shower or my position in it
    // pShower_ = nullptr;
    return *this;
  }
  
  Parton& Parton::operator=( const Parton &c)
  {
    JetScapeParticleBase::operator=(c);
    form_time_ = c.form_time_;
      Color_ = c.Color_;
      antiColor_ = c.antiColor_;
    set_edgeid ( c.edgeid() );
    pShower_ = c.shower();
    // set_edgeid( -1 ); // by default do NOT copy the shower or my position in it
    // pShower_ = nullptr;
    return *this;
  }

  void Parton::set_mean_form_time ()
  {
    mean_form_time_ = 2.0*e()/(t()+0.001);
  }
  
  void Parton::set_form_time(double form_time)
  {
    form_time_ = form_time;
  }
    
  void Parton::initialize_form_time()
  {
    form_time_ = -0.1;
  }

  double Parton::form_time()
  {
    return(form_time_);
  }
  
  const double Parton::mean_form_time()
  {
    return(mean_form_time_);
  }

  const double Parton::t()
  {
    /// \Todo: Fix 
    return ( PseudoJet::m2() ) ;
    // return (t_) ;
  }        

  void Parton::set_t(double t)
  {
    // This function has a very specific purpose and shouldn't normally be used
    // It scales down p! So catch people trying.
    if ( form_time()>= 0.0 ){
      throw std::runtime_error("Trying to set virtuality on a normal parton. You almost certainly don't want to do that. Please contact the developers if you do.");
    }
    
    //  Reset the momentum due to virtuality              
    double newPl = std::sqrt( e()*e() - t ) ;
    double velocityMod = std::sqrt(std::pow(jet_v_.comp(1),2) + std::pow(jet_v_.comp(2),2) + std::pow(jet_v_.comp(3),2));
    
    newPl = newPl/velocityMod;
    // double newP[4];
    // newP[0] = e();
    // for(int j=1;j<=3;j++) {
    //   newP[j] = newPl*jet_v_.comp(j);
    // }
    reset_momentum( newPl*jet_v_.comp(1), newPl*jet_v_.comp(2), newPl*jet_v_.comp(3), e() );
  } 

    void Parton::set_color(unsigned int col)
    {
        Color_ = col;
    }
    
    void Parton::set_anti_color(unsigned int acol)
    {
        antiColor_ = acol;
    }
    
    void Parton::set_max_color(unsigned int col)
    {
        MaxColor_ = col;
    }
    
    void Parton::set_min_color(unsigned int col)
    {
        MinColor_ = col;
    }
    
    void Parton::set_min_anti_color(unsigned int acol)
    {
        MinAntiColor_ = acol;
    }

  const int Parton::edgeid() const
  {
    return(edgeid_);
  }

  void Parton::set_edgeid( const int id )
  {
    edgeid_ = id;
  }

  void Parton::set_shower(const shared_ptr<PartonShower> pShower) {
    pShower_ = pShower;
  }
  
  const shared_ptr<PartonShower> Parton::shower() const{
    return pShower_;
  }

  std::vector<Parton> Parton::parents(){
    std::vector<Parton> ret;
    if ( !pShower_ ) return ret;
    node root = pShower_->GetEdgeAt(edgeid_).source();
    for ( node::in_edges_iterator parent = root.in_edges_begin(); parent != root.in_edges_end(); ++parent ){
      ret.push_back ( *pShower_->GetParton(*parent) );
    }
    return ret;
  }

    unsigned int Parton::color()
    {
        return (Color_);
    }
    
    unsigned int Parton::anti_color()
    {
        return (antiColor_);
    }
    
    unsigned int Parton::min_color()
    {
        return (MinColor_);
    }
    
    unsigned int Parton::min_anti_color()
    {
        return (MinAntiColor_);
    }

    unsigned int Parton::max_color()
    {
        return (MaxColor_);
    }
    
    // ---------------
    // Hadron specific
    // ---------------
    
    Hadron::Hadron (const Hadron& srh) :
    JetScapeParticleBase::JetScapeParticleBase (srh)
    {
        width_ = srh.width_ ;
    }
    
    Hadron::Hadron (int label, int id, int stat, const FourVector& p, const FourVector& x)  :
    JetScapeParticleBase::JetScapeParticleBase ( label,  id,  stat,  p, x)
    {
        set_decay_width(0.1);
    }
    
    
    Hadron::Hadron (int label, int id, int stat, double pt, double eta, double phi, double e, double* x)  :
    JetScapeParticleBase::JetScapeParticleBase ( label,  id,  stat,  pt, eta, phi, e, x)
    {
        set_decay_width(0.1);
        // cout << "========================== phieta Ctor called, returning : " << endl << *this << endl;
    }

    Hadron& Hadron::operator=( Hadron &c)
    {
        JetScapeParticleBase::operator=(c);
        width_ = c.width_;
        return *this;
    }
    
    Hadron& Hadron::operator=( const Hadron &c)
    {
        JetScapeParticleBase::operator=(c);
        width_ = c.width_;
        return *this;
    }

    
} /// end of namespace Jetscape
