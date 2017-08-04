// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// quick and dirty test class for Eloss modules ...
// can be used as a user template ...

#ifndef ELOSSMODULESTEST_H
#define ELOSSMODULESTEST_H

#include "JetEnergyLossModule.h"

// DEBUG ONLY
#include <fstream>
#include <sstream>

using namespace Jetscape;

class Matter : public JetEnergyLossModule<Matter> //, public std::enable_shared_from_this<Matter>
{  
 public:
  
  Matter();
  virtual ~Matter();

  void Init();
  //void Exec();
  //void DoEnergyLoss(double deltaT, double Q2, const vector<Parton>& pIn, vector<Parton>& pOut);
  void DoEnergyLoss(double deltaT, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);

  double generate_L(double form_time);
  double sudakov_Pgg(double g0, double g1, double loc_c, double E);
  double sud_val_GG(double h0, double h1, double h2, double loc_d, double E1);
  double sud_z_GG(double cg, double cg1, double loc_e , double l_fac, double E2);
  double P_z_gg_int(double cg, double cg1, double loc_e, double cg3, double l_fac , double E2);
  double sudakov_Pqg(double g0, double g1, double loc_c, double E);
  double sud_val_QG(double h0, double h1, double h2, double loc_d, double E1);
  double sud_z_QG(double cg, double cg1, double loc_e, double l_fac,double E2);
  double P_z_qg_int(double cg, double cg1, double loc_e, double cg3, double l_fac , double E2 );
  double sudakov_Pqq(double q0, double q1, double loc_c, double E);
  double sud_val_QQ(double h0, double h1, double h2, double loc_d, double E1);
  double sud_z_QQ(double cg, double cg1, double loc_e , double l_fac, double E2);
  double P_z_qq_int(double cg, double cg1, double loc_e, double cg3, double l_fac , double E2);
  //  void shower_vac( int line, int pid, double nu_in, double t0_in, double t_in, double kx, double ky, double loc, bool is_lead);
  double generate_vac_t(int p_id,double nu, double t0, double t, double loc_a, int isp);
  double  generate_vac_z(int p_id, double t0, double t, double loc_b, double nu, int is  );
  double alpha_s(double q2);
  double profile(double zeta);

  double generate_angle();
    
  double qhat,length;

  
 protected:
  uniform_real_distribution<double> ZeroOneDistribution;

  // DEBUG ONLY  
  ofstream* my_file_;

  
};


class Martini : public JetEnergyLossModule<Martini> //, public std::enable_shared_from_this<Martini>
{  
 public:
  
  Martini();
  virtual ~Martini();

  void Init();
  //void Exec();
  //void DoEnergyLoss(double deltaT, double Q2, const vector<Parton>& pIn, vector<Parton>& pOut);
  void DoEnergyLoss(double deltaT, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w) {};
  
 protected:
  uniform_real_distribution<double> ZeroOneDistribution;

  // DEBUG ONLY  
  ofstream* my_file_;

};

#endif

