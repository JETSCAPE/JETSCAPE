/** Provides JetScapeParticleBase and derived classes Parton, Hadron

    @group JetScape (modular/task) based framework
    @author Kolja Kauder
    @version Revision 0.1
    @date Aug 8, 2017
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

namespace Jetscape {


  // class Parton;
  // class Vertex;
  // class FourVector;



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
    double mass_            ; ///< rest mass of the particle \todo: only maintain PID, look up mass from PDG
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
