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

#ifndef MATTER_H
#define MATTER_H

#include "JetEnergyLossModule.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

class Matter : public JetEnergyLossModule<
                   Matter> //, public std::enable_shared_from_this<Matter>
{
public:
  Matter();
  virtual ~Matter();

  void Init();
  //void Exec();
  //void DoEnergyLoss(double deltaT, double Q2, const vector<Parton>& pIn, vector<Parton>& pOut);
  void DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton> &pIn,
                    vector<Parton> &pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);
  void Dump_pIn_info(int i, vector<Parton> &pIn);

  double generate_L(double form_time);
  double sudakov_Pgg(double g0, double g1, double loc_c, double E);
  double sud_val_GG(double h0, double h1, double h2, double loc_d, double E1);
  double sud_z_GG(double cg, double cg1, double loc_e, double l_fac, double E2);
  double P_z_gg_int(double cg, double cg1, double loc_e, double cg3,
                    double l_fac, double E2);
  double sudakov_Pqg(double g0, double g1, double loc_c, double E);
  double sud_val_QG(double h0, double h1, double h2, double loc_d, double E1);
  double sud_z_QG(double cg, double cg1, double loc_e, double l_fac, double E2);
  double P_z_qg_int(double cg, double cg1, double loc_e, double cg3,
                    double l_fac, double E2);
  double sudakov_Pqg_w_M(double M, double g0, double g1, double loc_c,
                         double E);
  double sud_val_QG_w_M(double M, double h0, double h1, double h2, double loc_d,
                        double E1);
  double sud_z_QG_w_M(double M, double cg, double cg1, double loc_e,
                      double l_fac, double E2);
  double P_z_qg_int_w_M(double M, double cg, double cg1, double loc_e,
                        double cg3, double l_fac, double E2);
  double sudakov_Pqq(double q0, double q1, double loc_c, double E);

  double sud_val_QQ(double h0, double h1, double h2, double loc_d, double E1);
  double sud_z_QQ(double cg, double cg1, double loc_e, double l_fac, double E2);
  double P_z_qq_int(double cg, double cg1, double loc_e, double cg3,
                    double l_fac, double E2);
  double P_z_qp_int(double cg, double cg1, double loc_e, double cg3,
                    double l_fac, double E2);
  double sud_z_QP(double cg, double cg1, double loc_e, double l_fac, double E2);
  double sud_val_QP(double h0, double h1, double h2, double loc_d, double E1);
  double sudakov_Pqp(double g0, double g1, double loc_c, double E);

  double sudakov_Pqq_w_M_vac_only(double M, double q0, double q1, double loc_c,
                                  double E);
  double sud_val_QQ_w_M_vac_only(double M, double h0, double h1, double h2,
                                 double loc_d, double E1);
  double sud_z_QQ_w_M_vac_only(double M, double cg, double cg1, double loc_e,
                               double l_fac, double E2);
  double P_z_qq_int_w_M_vac_only(double M, double cg, double cg1, double loc_e,
                                 double cg3, double l_fac, double E2);

  //  void shower_vac( int line, int pid, double nu_in, double t0_in, double t_in, double kx, double ky, double loc, bool is_lead);
  double generate_vac_t(int p_id, double nu, double t0, double t, double loc_a,
                        int isp);
  double generate_vac_t_w_M(int p_id, double M, double nu, double t0, double t,
                            double loc_a, int is);
  double generate_vac_z(int p_id, double t0, double t, double loc_b, double nu,
                        int is);
  double generate_vac_z_w_M(int p_id, double M, double t0, double t,
                            double loc_b, double nu, int is);
  double alpha_s(double q2);
  double profile(double zeta);

  double generate_angle();
  double generate_kt(double local_qhat, double dzeta);

  double qhat = 0.0;
  double ehat = 0.0;
  double e2hat = 0.0;
  double length = 0.0;

  unsigned int MaxColor = 0;

  //SC: for interface with hydro
  double fillQhatTab(double y);
  double fncQhat(double zeta);
  double fncAvrQhat(double zeta, double tau);

  bool matter_on, in_vac, brick_med, recoil_on, broadening_on;
  double hydro_Tc, qhat0, alphas, brick_length, vir_factor;
  double initR0, initRx, initRy, initRz, initVx, initVy, initVz, initRdotV,
      initVdotV, initEner;
  double Q00, Q0, T0;

  static const int dimQhatTab = 151;
  double qhatTab1D[dimQhatTab] = {0.0};
  double qhatTab2D[dimQhatTab][dimQhatTab] = {{0.0}};

  double tStart = 0.6;
  int iEvent;
  bool debug_flag = 0;
  long NUM1;

  // Variables for HQ 2->2
  static const int N_p1 = 500;
  static const int N_T = 60;
  static const int N_e2 = 75;
  static double distFncB[N_T][N_p1][N_e2], distFncF[N_T][N_p1][N_e2],
      distMaxB[N_T][N_p1][N_e2], distMaxF[N_T][N_p1][N_e2];
  static double distFncBM[N_T][N_p1], distFncFM[N_T][N_p1];
  double min_p1 = 0.0;
  double max_p1 = 1000.0;
  double bin_p1 = (max_p1 - min_p1) / N_p1;
  double min_T = 0.1;
  double max_T = 0.7;
  double bin_T = (max_T - min_T) / N_T;
  double min_e2 = 0.0;
  double max_e2 = 15.0;
  double bin_e2 = (max_e2 - min_e2) / N_e2;

  static double RHQ[60][20];    //total scattering rate for heavy quark
  static double RHQ11[60][20];  //Qq->Qq
  static double RHQ12[60][20];  //Qg->Qg
  static double qhatHQ[60][20]; //qhat of heavy quark

  // flag to make sure initialize only once
  static bool flag_init;

  //SC: for elastic scattering
  void flavor(int &CT, int &KATT0, int &KATT2, int &KATT3,
              unsigned int &max_color, unsigned int &color0,
              unsigned int &anti_color0, unsigned int &color2,
              unsigned int &anti_color2, unsigned int &color3,
              unsigned int &anti_color3);
  void colljet22(int CT, double temp, double qhat0ud, double v0[4],
                 double p0[4], double p2[4], double p3[4], double p4[4],
                 double &qt);
  void trans(double v[4], double p[4]);
  void transback(double v[4], double p[4]);
  void rotate(double px, double py, double pz, double pr[4], int icc);
  float ran0(long *idum);
  double solve_alphas(double var_qhat, double var_ener, double var_temp);
  double fnc0_alphas(double var_alphas, double var_qhat, double var_ener,
                     double var_temp);
  double fnc0_derivative_alphas(double var_alphas, double var_qhat,
                                double var_ener, double var_temp);

  void read_tables();
  double Mgc2gc(double s, double t, double M);
  double Mqc2qc(double s, double t, double M);
  void collHQ22(int CT, double temp, double qhat0ud, double v0[4], double p0[4],
                double p2[4], double p3[4], double p4[4], double &qt);

protected:
  uniform_real_distribution<double> ZeroOneDistribution;

private:
  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<Matter> reg;
};

#endif // MATTER_H
