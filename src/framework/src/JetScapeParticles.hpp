/** Provides JetScapeParticleBase and derived classes Parton, Hadron

    @group JetScape (modular/task) based framework
    @author Kolja Kauder
    @version Revision 0.1
    @date Aug 8, 2017

    \class Jetscape::JetScapeParticleBase
    * A JetScapeParticleBase is a FastJet PseudoJet with additional information
    *  - PID (from PDG) and rest mass (these should eventually be coupled and only PID kept track of internally)
    *  - A location (creation point) 4-vector
    *  - a label and a status 
    *  - currently additional information that should be moved to derived classes or only used as UserInfo
    * 
    * You can in principle use the base class directly, but it's recommended to use the derived classes
    * Parton and/or (todo) Hadron, Lepton, ...
    * 
    * \warning
    * PseudoJet doesn't have a concept of rest mass, any mass related functions literally assume 
    * \f$M^2 = E^2 - p*p.\f$
    * Especially in the case of off-shell partons, the correct interpretation is 
    * \f$"M^2" = E^2 - p^2 == M0^2 + Q^2 == M0^2 + t\f$
    * Therefore, there is the dangerous possibility functions like mass(), reset_PtYPhiM(...) can
    * be used by unwitting users and display unexpected behavior. 
    * For protection against this possibility, "mass"-related functions in PseudoJet are
    * overwritten to throw an error.
    * This CAN be circumvented, like
    \code{.cpp}
    JetScapeParticleBase j(...);
    PseudoJet* pj = (PseudoJet*) &j;
    cout << j->mass() << endl;
    \endcode
    * Since PseudoJet methods are not virtual there is no way to prohibit this, 
    * but it can be considered a good thing for experienced users who know what they're doing.
    * 
    * Future considerations: 
    *   - We should consider
    *     making a Pythia8 installation mandatory; with pythia guaranteed to be present,
    *     the rest mass lookup could be automatically done using PDG data.
    *   - If ROOT were a mandatory part, TLorentzVector would be a good replacement for the homebrewed FourVector
*/

#ifndef JetScapeParticles_hpp
#define JetScapeParticles_hpp

#include <stdio.h>
#include <math.h>
#include "constants.h"
#include "four_vector.hpp"
#include "fjcore.hh"

#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>

/** 
*/

namespace Jetscape {

  /**************************************************************************************************/
  //  BASE CLASS
  /*************************************************************************************************/

  class JetScapeParticleBase : public fjcore::PseudoJet
  {

  public:

    JetScapeParticleBase() : PseudoJet() {};
    JetScapeParticleBase (int label, int id, int stat, double p[4], double x[4]);
    JetScapeParticleBase (int label, int id, int stat, double pt, double eta, double phi, double e);
    JetScapeParticleBase (int label, int id, int stat, double pt, double eta, double phi, double e, double x[4]);
    JetScapeParticleBase (const JetScapeParticleBase& srp);
	  
    virtual ~JetScapeParticleBase();
  
    void clear();    

    // Setters
    void set_label(int label);
    void set_id(int id);  
    void set_stat(int stat);

    // from PseudoJet ... (well keep for now);
    void set_mass(double mass_input);    
    void set_p(double p[4]);
    void set_x(double x[4]); 
    
    void set_mean_form_time();
    void set_form_time(double form_time);    
    void set_t(double t); ///< virtuality of particle
    
    void init_jet_v();
    void set_jet_v(double v[4]);
    
    //  Getters
    
    double form_time();
    const int pid();    
    const int pstat();    
    const int plabel();
    const double e();
    const double pt();
    const double time();
    
    FourVector &p_in();  
    FourVector &x_in();
    FourVector &jet_v();
  
    const double mass();
    const double mean_form_time();
    const double p(int i); 
    double pl();    
    const double nu();    
    const double t();    
    const double t_max();

    JetScapeParticleBase& operator=(JetScapeParticleBase &c);
    JetScapeParticleBase& operator=(const JetScapeParticleBase &c);

  protected:
  
    int pid_                ; ///< particle id ()
    int pstat_              ; ///< status of particle
    int plabel_             ; ///< the line number in the event record
    double mass_            ; ///< rest mass of the particle \todo Only maintain PID, look up mass from PDG
    double t_               ; ///< The virtuality, and not the time!
    double mean_form_time_  ; ///< Mean formation time
    double form_time_       ; ///< event by event formation time
    
    FourVector p_in_; ///< Internal version of p. Clashes with PseudoJet! REPLACE
    FourVector x_in_; ///< position of particle
    FourVector jet_v_; ///< jet four vector, without gamma factor (so not really a four vector)

    // helpers
    void initialize_form_time();
  };
  // END BASE CLASS


  // Declared outside the class
  ostream &operator<<( ostream &output, JetScapeParticleBase & p );

  /**************************************************************************************************/
  //  PARTON CLASS
  /*************************************************************************************************/
  class Parton : public JetScapeParticleBase{
    using JetScapeParticleBase::JetScapeParticleBase;
  };


};  /// end of namespace Jetscape

#endif /* JetScapeParticles */
