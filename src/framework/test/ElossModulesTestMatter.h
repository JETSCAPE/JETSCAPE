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

using namespace Jetscape;

class Matter : public JetEnergyLossModule<Matter> //, public std::enable_shared_from_this<Matter>
{  
 public:
  
  Matter();
  virtual ~Matter();

  void Init();
  //void Exec();
  //void DoEnergyLoss(double deltaT, double Q2, const vector<Parton>& pIn, vector<Parton>& pOut);
  void DoEnergyLoss(double deltaT,double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
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

  //SC: for interface with hydro
  double fillQhatTab();
  double fncQhat(double zeta);
  double fncAvrQhat(double zeta, double tau);

  bool in_vac,brick_med;
  double hydro_Tc,qhat0,alphas,brick_length,vir_factor;
  double initR0,initRx,initRy,initRz,initVx,initVy,initVz,initRdotV,initEner;
  double Q00,Q0,T0;

  static const int dimQhatTab=151;
  double qhatTab1D[dimQhatTab]={0.0};
  double qhatTab2D[dimQhatTab][dimQhatTab]={{0.0}};

  long  NUM1;                  

  //SC: for elastic scattering
  void flavor(int &CT,int &KATT0,int &KATT2,int &KATT3, unsigned int &max_color, unsigned int &color0, unsigned int &anti_color0, unsigned int &color2, unsigned int &anti_color2, unsigned int &color3, unsigned int &anti_color3);
  void colljet22(int CT,double temp,double qhat0ud,double v0[4],double p0[4],double p2[4],double p3[4],double p4[4],double &qt);
  void trans(double v[4],double p[4]);
  void transback(double v[4],double p[4]);
  void rotate(double px,double py,double pz,double pr[4],int icc);
  float ran0(long *idum);
  double solve_alphas(double var_qhat, double var_ener, double var_temp);
  double fnc0_alphas(double var_alphas, double var_qhat, double var_ener, double var_temp);
  double fnc0_derivative_alphas(double var_alphas, double var_qhat, double var_ener, double var_temp);

  
 protected:
  uniform_real_distribution<double> ZeroOneDistribution;
  
};


#endif

