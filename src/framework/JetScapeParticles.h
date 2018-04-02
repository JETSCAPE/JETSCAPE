/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/amajumder/JETSCAPE-COMP/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

/** Provides JetScapeParticleBase and derived classes Parton, Hadron

    @group JetScape (modular/task) based framework
    @author Kolja Kauder
    @version Revision 0.1
    @date Aug 8, 2017

    \class Jetscape::JetScapeParticleBase
    * A JetScapeParticleBase derives PRIVARTELY from FastJet PseudoJet and has additional information:
    *  - PID (from PDG) and rest mass (these should eventually be coupled and only PID kept track of internally)
    *  - A location (creation point) 4-vector
    *  - a label and a status 
    *  - currently additional information that should be moved to derived classes or only used as UserInfo
    * 
    * The design choice of protected inheritance is due to a disconnect between available packages.
    * The overwhelming majority of the theory community expects the 0 component to be time/energy, 
    * whereas FastJet (and others, like ROOT) prefer t,e to be the fourth component.
    * Private inheritance means we can inherit and make accessible safe methods (with C++11 using),
    * while protecting users from unsafe (explicit component access) ones.
    * Note that this is only necessary because otherwise it's impossible to disallow 
    * constructors and getters that explicitly assume indices!
    * IF we could get rid of those or change to the fastjet convention,
    * we could derive publicly and get a true Is_A relationship. Alas.
    * 
    * You can in principle use the base class directly, but it's recommended to use the derived classes
    * Parton and/or (todo) Hadron, Lepton, ...
    * 
    * \warning
    * PseudoJet doesn't have a concept of rest mass, any mass related functions literally assume 
    * \f$M^2 = E^2 - p*p.\f$
    * Especially in the case of off-shell partons, the correct interpretation is 
    * \f$"M^2" = E^2 - p^2 == M0^2 + Q^2 == M0^2 + t\f$
    * Therefore, there is the dangerous possibility that functions like mass(), reset_PtYPhiM(...) can
    * be used by unwitting users and display unexpected behavior. 
    * For protection against this possibility, "mass"-related functions in PseudoJet are
    * not made available
    * 
    * Future considerations: 
    *   - We should consider
    *     making a Pythia8 installation mandatory; with pythia guaranteed to be present,
    *     the rest mass lookup could be automatically done using PDG data.
    *   - If ROOT were a mandatory part, TLorentzVector would be a good replacement for the homebrewed FourVector
    *   - If HepMc were a mandatory part, HepMc::FourVector would be a good replacement for the homebrewed FourVector
    */

#ifndef JETSCAPEPARTICLES_H
#define JETSCAPEPARTICLES_H

#include <stdio.h>
#include <math.h>
#include "GTL/graph.h"
#include "JetScapeConstants.h"
#include "FourVector.h"
#include "fjcore.hh"
#include "JetScapeLogger.h"
#include "PartonShower.h"

#include "Pythia8/Pythia.h"

#include <vector>
#include <memory>
#include <iostream>
#include <sstream>
#include <iomanip>

using std::ostream;
using std::weak_ptr;

namespace Jetscape {

  class PartonShower;
  
  /**************************************************************************************************/
  //  BASE CLASS
  /*************************************************************************************************/

  class JetScapeParticleBase : protected fjcore::PseudoJet
  {
    friend class fjcore::PseudoJet;

    // unsafe
    // using fjcore::PseudoJet::PseudoJet;
    // using fjcore::PseudoJet::operator() (int i) const ; 
    // inline double operator [] (int i) const { return (*this)(i); }; // this too

  public:
    // Disallow reset, it assumes different logic
    // using fjcore::PseudoJet::reset(double px, double py, double pz, double E);
    // using fjcore::PseudoJet::reset(const PseudoJet & psjet) {
    //   inline void reset_momentum(double px, double py, double pz, double E);
    //   inline void reset_momentum(const PseudoJet & pj);

    inline void reset_momentum( const double px, const double py, const double pz, const double e ){
      fjcore::PseudoJet::reset_momentum ( px, py, pz, e );
    }

    inline void reset_momentum( const FourVector& p ){
      fjcore::PseudoJet::reset_momentum ( p.x(), p.y(), p.z(), p.t() );
    }
    
    // Disallow the valarray return.
    // Can replace/and or provide a double[4] version with the right assumptions
    // std::valarray<double> four_mom() const;
    // enum { X=0, Y=1, Z=2, T=3, NUM_COORDINATES=4, SIZE=NUM_COORDINATES };

    // This _should_ work, but not taking chances for now
    // PseudoJet & boost(const PseudoJet & prest);
    // PseudoJet & unboost(const PseudoJet & prest);

    // Replace with appropriate functions. m is a tricky one. 
    // inline double  m2() const {return (_E+_pz)*(_E-_pz)-_kt2;}    
    // inline double  m() const;    
    // inline double mperp2() const {return (_E+_pz)*(_E-_pz);}
    // inline double mperp() const {return sqrt(std::abs(mperp2()));}
    // inline double mt2() const {return (_E+_pz)*(_E-_pz);}
    // inline double mt() const {return sqrt(std::abs(mperp2()));}

    // Disallow functions containing M
    // inline void reset_PtYPhiM(...);
    // void reset_momentum_PtYPhiM(double pt, double y, double phi, double m=0.0);

    // // void set_cached_rap_phi(double rap, double phi);

    /// No implicit cast to PseudoJet is allowed, provide a conversion
    fjcore::PseudoJet GetPseudoJet() const{
      return PseudoJet ( *this );
    }
    
    // import safe functions
    using fjcore::PseudoJet::px;
    using fjcore::PseudoJet::py;
    using fjcore::PseudoJet::pz;
    using fjcore::PseudoJet::e;
    using fjcore::PseudoJet::E;

    using fjcore::PseudoJet::phi;
    using fjcore::PseudoJet::phi_std;
    using fjcore::PseudoJet::phi_02pi;
    using fjcore::PseudoJet::rap;
    using fjcore::PseudoJet::rapidity;
    using fjcore::PseudoJet::pseudorapidity;
    using fjcore::PseudoJet::eta;
    using fjcore::PseudoJet::pt2;
    using fjcore::PseudoJet::pt;
    using fjcore::PseudoJet::perp2;
    using fjcore::PseudoJet::perp;
    using fjcore::PseudoJet::kt2;

    using fjcore::PseudoJet::modp2;
    using fjcore::PseudoJet::modp;
    using fjcore::PseudoJet::Et;
    using fjcore::PseudoJet::Et2;

    using fjcore::PseudoJet::kt_distance;
    using fjcore::PseudoJet::plain_distance;
    using fjcore::PseudoJet::squared_distance;
    using fjcore::PseudoJet::delta_R;
    using fjcore::PseudoJet::delta_phi_to;
    using fjcore::PseudoJet::beam_distance;

    using fjcore::PseudoJet::operator*=;
    using fjcore::PseudoJet::operator/=;
    using fjcore::PseudoJet::operator+=;
    using fjcore::PseudoJet::operator-=;

    using fjcore::PseudoJet::user_index;
    using fjcore::PseudoJet::set_user_index;
    using fjcore::PseudoJet::UserInfoBase;
    using fjcore::PseudoJet::InexistentUserInfo;

    using fjcore::PseudoJet::user_info;
    using fjcore::PseudoJet::set_user_info;
    using fjcore::PseudoJet::has_user_info;
    using fjcore::PseudoJet::user_info_ptr;
    using fjcore::PseudoJet::user_info_shared_ptr;

    using fjcore::PseudoJet::description;
    // In principle, these might be okay, but ClusterSequences should
    // be made after explicitly transforming to a proper PseudoJet
    // using fjcore::PseudoJet::has_associated_cluster_sequence;
    // using fjcore::PseudoJet::has_associated_cs;
    // using fjcore::PseudoJet::has_valid_cluster_sequence;
    // using fjcore::PseudoJet::has_valid_cs;
    // using fjcore::PseudoJet::associated_cluster_sequence;
    // using fjcore::PseudoJet::associated_cs;
    // using fjcore::PseudoJet::validated_cluster_sequence;
    // using fjcore::PseudoJet::validated_cs;
    // using fjcore::PseudoJet::set_structure_shared_ptr;
    // using fjcore::PseudoJet::has_structure;
    // using fjcore::PseudoJet::structure_ptr;
    // using fjcore::PseudoJet::structure_non_const_ptr;
    // using fjcore::PseudoJet::validated_structure_ptr;
    // using fjcore::PseudoJet::structure_shared_ptr;
    // ... more
    
  public:

    JetScapeParticleBase() : PseudoJet() {};
    JetScapeParticleBase (int label, int id, int stat, const FourVector& p, const FourVector& x);
    JetScapeParticleBase (int label, int id, int stat, double pt, double eta, double phi, double e, double* x=0);
    JetScapeParticleBase(int label, int id, int stat,
                         const FourVector& p, const FourVector& x,
                         double mass);
    JetScapeParticleBase (const JetScapeParticleBase& srp);
	  
    virtual ~JetScapeParticleBase();
  
    void clear();    

    // Setters
    void set_label(int label);
    void set_id(int id);  
    void set_stat(int stat);
    void set_x(double x[4]); 
    
    void init_jet_v();
    void set_jet_v(double v[4]);
    void set_jet_v(FourVector j);
    
    //  Getters
    
    const int pid() const;
    const int pstat() const;
    const int plabel() const;
    // const double e();
    // const double pt();
    const double time() const;
    
    std::vector<JetScapeParticleBase> parents();

    const FourVector p_in() const;
    const FourVector &x_in() const;
    const FourVector &jet_v() const;

    const double restmass();
    const double p(int i);
    double pl();
    const double nu();
    const double t_max();
    
    virtual JetScapeParticleBase& operator=(JetScapeParticleBase &c);
    virtual JetScapeParticleBase& operator=(const JetScapeParticleBase &c);

  protected:
  
    void set_restmass(double mass_input); ///< shouldn't be called from the outside, needs to be consistent with PID

    int pid_                ; ///< particle id
    int pstat_              ; ///< status of particle
    int plabel_             ; ///< the line number in the event record
    double mass_            ; ///< rest mass of the particle \todo Only maintain PID, look up mass from PDG

    FourVector x_in_; ///< position of particle
    FourVector jet_v_; ///< jet four vector, without gamma factor (so not really a four vector)

    // give it a static pythia to look up particle properties
    // Be a bit careful with it!
    // Init is never called, and this object is not configured. All it can do is look up
    // in its original Data table
    static Pythia8::Pythia InternalHelperPythia;  

  };
  // END BASE CLASS


  // Declared outside the class
  ostream &operator<<( ostream &output, JetScapeParticleBase & p );

  /**************************************************************************************************/
  //  PARTON CLASS
  /*************************************************************************************************/
  class Parton : public JetScapeParticleBase{
    
  public :
    virtual void set_mean_form_time();
    virtual void set_form_time(double form_time);    

    virtual double form_time();
    virtual const double mean_form_time();
    virtual void set_color(unsigned int col); ///< sets the color of the parton
    virtual void set_anti_color(unsigned int acol); ///< sets anti-color of the parton
    virtual void set_max_color(unsigned int col); ///< sets the color of the parton
    virtual void set_min_color(unsigned int col); ///< sets the color of the parton
    virtual void set_min_anti_color(unsigned int acol); ///< sets anti-color of the parton

    Parton (int label, int id, int stat, const FourVector& p, const FourVector& x);
    Parton (int label, int id, int stat, double pt, double eta, double phi, double e, double* x=0);
    Parton (const Parton& srp);
    
    Parton& operator=( Parton &c);
    Parton& operator=( const Parton &c);
    
    const double t();
    void set_t(double t); ///< virtuality of particle, \WARNING: rescales the spatial component
    unsigned int color(); ///< returns the color of the parton
    unsigned int anti_color(); ///< returns the anti-color of the parton
    unsigned int max_color();
    unsigned int min_color();
    unsigned int min_anti_color();

    
    const int edgeid() const;
    void set_edgeid( const int id);

    void set_shower(const shared_ptr<PartonShower> pShower);
    void set_shower(const weak_ptr<PartonShower> pShower);
    const weak_ptr<PartonShower> shower() const;

    std::vector<Parton> parents();
    
  protected :
    double mean_form_time_  ; ///< Mean formation time
    double form_time_       ; ///< event by event formation time
    unsigned int Color_      ; ///< Large Nc color of parton
    unsigned int antiColor_  ;///< Large Nc anti-color of parton
    unsigned int MaxColor_    ; ///< the running maximum color
    unsigned int MinColor_    ; ///< color of the parent
    unsigned int MinAntiColor_; ///< anti-color of the parent

    weak_ptr<PartonShower> pShower_; ///< shower that this parton belongs to
    int edgeid_             ; ///< Position in the shower graph    

    // helpers
    void initialize_form_time();
    void CheckAcceptability ( int id ); ///< restrict to a few pids only.
    
  };

  class Hadron : public JetScapeParticleBase
  {
    public:
        
    Hadron (int label, int id, int stat, const FourVector& p, const FourVector& x);
    Hadron (int label, int id, int stat, double pt, double eta, double phi, double e, double* x=0);
    Hadron(int label, int id, int stat,
           const FourVector& p, const FourVector& x, double mass);
    Hadron (const Hadron& srh);
        
    Hadron& operator=( Hadron &c);
    Hadron& operator=( const Hadron &c);

    void set_decay_width(double width)
    {
      width_ = width;
    }
        
    double decay_width()
    {
      return(width_);
    }

    /// Hadron may be used to handle electrons, gammas, ... as well
    /// In addition, not all generated ids may be in the database
    /// Currently, we add these manually. Could also reject outright.
    bool CheckOrForceHadron( const int id, const double mass=0 );
      
  protected:
    double width_;
        
  };
    
    
};  /// end of namespace Jetscape

#endif // JETSCAPEPARTICLES_H
