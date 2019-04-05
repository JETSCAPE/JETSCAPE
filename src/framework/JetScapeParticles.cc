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

// Provides JetScapeParticleBase and derived classes Parton, Hadron

#include <iostream>
#include <complex>
#include <fstream>
#include <cmath>
#include <assert.h>
#include "JetScapeLogger.h"
#include "JetScapeParticles.h"
#include "JetScapeConstants.h"

namespace Jetscape {

  // Initialize static MakeUniqueHelper.here
  Pythia8::Pythia JetScapeParticleBase::InternalHelperPythia ("IntentionallyEmpty",false);

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
    jet_v_ = srp.jet_v_;
    x_in_ = srp.x_in_;
  }

  JetScapeParticleBase::JetScapeParticleBase (int label, int id, int stat, double pt, double eta, double phi, double e, double* x)
  {
    set_label(label);
    set_id(id);
    init_jet_v();
    
    assert ( InternalHelperPythia.particleData.isParticle(id) );
    set_restmass( InternalHelperPythia.particleData.m0( id ) );
    
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
    
  }

    
  JetScapeParticleBase::JetScapeParticleBase (int label, int id, int stat, const FourVector& p, const FourVector& x)
  {
    
    set_label(label);
    set_id(id);    
    init_jet_v();
    
    assert ( InternalHelperPythia.particleData.isParticle(id) );
    set_restmass( InternalHelperPythia.particleData.m0( id ) );

    reset_momentum(p);
    x_in_=x;
    set_stat(stat);

  }
  
  JetScapeParticleBase::JetScapeParticleBase(
                    int label, int id, int stat,
                    const FourVector& p, const FourVector& x, double mass) {
    set_label(label);
    set_id(id);    
    init_jet_v();
    
    reset_momentum(p);
    x_in_=x;
    set_stat(stat);

  }


  void JetScapeParticleBase::clear()
  {
    plabel_ = 0;
    pid_ = 0;
    pstat_ = 0;
    mass_=-1;
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
  const int JetScapeParticleBase::pid() const
  {
    return(pid_);
  }
    
  const int JetScapeParticleBase::pstat() const
  {
    return(pstat_);
  }
    
  const int JetScapeParticleBase::plabel() const
  {
    return(plabel_);
  }

      
  const double JetScapeParticleBase::time() const
  {
    return(x_in_.t());
  }
      
  const FourVector JetScapeParticleBase::p_in() const
  {
    return( FourVector(px(), py(), pz(), e() ));
  }

  const FourVector &JetScapeParticleBase::x_in() const
  {
    return(x_in_);
  }
    
  const FourVector &JetScapeParticleBase::jet_v() const 
  {
    return(jet_v_);
  }
  
  const double JetScapeParticleBase::restmass()
  {
    return(mass_);
  }

  const double JetScapeParticleBase::p(int i) {
    /// Deprecated. Prefer explicit component access
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
    set_shower ( srp.shower() );

    // set_edgeid( -1 ); // by default do NOT copy the shower or my position in it
    // pShower_ = nullptr;
  }

  Parton::Parton (int label, int id, int stat, const FourVector& p, const FourVector& x)  :
    JetScapeParticleBase::JetScapeParticleBase ( label,  id,  stat,  p, x)
  {
    CheckAcceptability( id );
    assert ( InternalHelperPythia.particleData.isParton(id) );
    initialize_form_time();
    set_color(0);
    set_anti_color(0);
    set_min_color(0);
    set_min_anti_color(0);
    set_max_color(0);
    set_edgeid( -1 );
    set_shower ( 0 );
    
    // cout << "========================== std Ctor called, returning : " << endl << *this << endl;
  }
  

  Parton::Parton (int label, int id, int stat, double pt, double eta, double phi, double e, double* x)  :
    JetScapeParticleBase::JetScapeParticleBase ( label,  id,  stat,  pt, eta, phi, e, x){
    CheckAcceptability ( id );
    assert ( InternalHelperPythia.particleData.isParton(id) );
    initialize_form_time();
    set_color(0);
    set_anti_color(0);
    set_min_color(0);
    set_min_anti_color(0);
    set_max_color(0);
    set_edgeid( -1 );
    set_shower ( 0 );

    // cout << "========================== phieta Ctor called, returning : " << endl << *this << endl;
  }
  
  void Parton::CheckAcceptability ( int id ){
    switch (id) {
    case 1:  //down quark
    case -1: // anti-down quark
      break;            
    case 2:  // up quark
    case -2:  // anti-up quark
      break;
    case 3:   // strange quark
    case -3:  // anti-strange quark
      break;
    case 4:   // charm quark
    case -4:  // anti-charm quark
      break;            
    case 5:   // bottom quark
    case -5:  // anti-bottom quark
      break;            
    case 21: // gluon
      break;      
    default:
      JSWARN << " error in id = " << id;
      throw std::runtime_error ( "pid not accepted for Parton");
      break;
    }
  }


  Parton& Parton::operator=( Parton &c)
  {
    JetScapeParticleBase::operator=(c);
    form_time_ = c.form_time_;
    Color_ = c.Color_;
    antiColor_ = c.antiColor_;
    set_edgeid ( c.edgeid() );
    set_shower ( c.shower() );

    return *this;
  }
  
  Parton& Parton::operator=( const Parton &c)
  {
    JetScapeParticleBase::operator=(c);
    form_time_ = c.form_time_;
    Color_ = c.Color_;
    antiColor_ = c.antiColor_;
    set_edgeid ( c.edgeid() );
    set_shower ( c.shower() );

    return *this;
  }

  void Parton::set_mean_form_time ()
  {
    mean_form_time_ = 2.0*e()/(t()+0.001)/fmToGeVinv;
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

    void Parton::reset_p(double px, double py, double pz )
    {
        reset_momentum(px,py,pz,e());
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
    if ( pShower!= nullptr )
      pShower_ = pShower;
  }

  void Parton::set_shower(const weak_ptr<PartonShower> pShower) {
      pShower_ = pShower;
  }
  
  const weak_ptr<PartonShower> Parton::shower() const{
    return pShower_;
  }

  std::vector<Parton> Parton::parents(){
    std::vector<Parton> ret;
    auto shower = pShower_.lock();
    if ( !shower ) return ret;
    node root = shower->GetEdgeAt(edgeid_).source();
    for ( node::in_edges_iterator parent = root.in_edges_begin(); parent != root.in_edges_end(); ++parent ){
      ret.push_back ( *shower->GetParton(*parent) );
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
    assert ( CheckOrForceHadron( id ) );
    // assert ( InternalHelperPythia.particleData.isHadron(id) );
    set_decay_width(0.1);
  }
    
    
  Hadron::Hadron (int label, int id, int stat, double pt, double eta, double phi, double e, double* x)  :
    JetScapeParticleBase::JetScapeParticleBase ( label,  id,  stat,  pt, eta, phi, e, x)
  {
    assert ( CheckOrForceHadron( id ) );
    // assert ( InternalHelperPythia.particleData.isHadron(id) );
    set_decay_width(0.1);
    // cout << "========================== phieta Ctor called, returning : " << endl << *this << endl;
  }
  
  Hadron::Hadron(int label, int id, int stat,
                 const FourVector& p, const FourVector& x, double mass):
    JetScapeParticleBase::JetScapeParticleBase(label, id, stat, p, x, mass) {
    assert ( CheckOrForceHadron( id, mass ) );
    set_restmass(mass);
  }

  bool Hadron::CheckOrForceHadron( const int id, const double mass ){
    bool status = InternalHelperPythia.particleData.isHadron(id);
    if (status ) return true;
    
    // If it's not recognized as a hadron, still allow some (or all)
    // particles. Particularly leptons and gammas are the point here.
    // TODO: Handle non-partonic non-hadrons more gracefully

    // -- Add unknown particles
    if  ( !InternalHelperPythia.particleData.isParticle(id) ) { // avoid doing it over and over
      JSWARN << "id = " << id << " is not recognized as a hadron! "
	   << "Add it as a new type of particle.";
      InternalHelperPythia.particleData.addParticle( id, " ", 0, 0, 0, mass, 0.1);
    }

    // -- now all that's left is known non-hadrons. We'll just accept those.
    return true;      
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
