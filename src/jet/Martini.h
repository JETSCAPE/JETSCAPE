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

#ifndef MARTINI_H
#define MARTINI_H

#include <fstream>
#include <math.h>
#include "JetEnergyLossModule.h"
#include "JetScapeConstants.h"

using namespace Jetscape;

//Basic.h//
struct RateRadiative
{
  double qqg;
  double ggg;
  double gqq;
  double qqgamma;
};

struct RateElastic
{
  double qq;
  double gq;
  double qg;
  double gg;
};

struct RateConversion
{
  double qg;
  double gq;
  double qgamma;
};

class Martini : public JetEnergyLossModule<Martini> //, public std::enable_shared_from_this<Martini>
{  
 private:

  // AMY rates are calculated in p/T > AMYpCut
  static constexpr double AMYpCut = 4.01;

  double alpha_s;
  double alpha_em;
  double g;
  double pcut;        // below this scale, no further Eloss

  //Import.h//
  static const int NP = 230;
  static const int NK = 381;

  static const int Nalphas = 11;
  static const int Nomega = 120;
  static const int Nq = 60;

  static constexpr double omegaStep = 0.2;
  static constexpr double qStep = 0.2;
  static constexpr double alphaMin = 0.15;
  static constexpr double alphaStep = 0.03;

  typedef struct
  {
    double ddf;
    double dda;
    double dcf;
    double dca;
    int    include_gluons;
    int Nc;
    int Nf;
    int BetheHeitler;
    int BDMPS;
    double alpha;
    double delta_x;
    double dx_max;
    double dp;
    double p_min;
    double p_max;
    long   n_p;
    long   n_pmin;
    double k_min;
    double k_max;
    long   n_k;
    long   n_kmin;
  } Gamma_info;

  // Structure to store information about splitting functions
  typedef struct
  {
    double qqg[NP][NK];
    double gqq[NP][NK];
    double ggg[NP][NK];
    double qqgamma[NP][NK];

    double tau_qqg[NP][NK];
    double tau_gqq[NP][NK];
    double tau_ggg[NP][NK];
    double tau_qqgamma[NP][NK];
  } dGammas;

  Gamma_info dat;
  dGammas    Gam;

  vector<double> *dGamma_qq;
  vector<double> *dGamma_qg;
  vector<double> *dGamma_qq_q;
  vector<double> *dGamma_qg_q;

 public:
  
  Martini();
  virtual ~Martini();

  //main//
  void Init();
  void DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  int DetermineProcess(double p, double T, double deltaT, int id);
  void WriteTask(weak_ptr<JetScapeWriter> w) {};
  
  //Radiative.h//
  RateRadiative getRateRadTotal(double p, double T);
  RateRadiative getRateRadPos(double u, double T);
  RateRadiative getRateRadNeg(double u, double T);

  double getNewMomentumRad(double p, double T, int process);
  double area(double y, double u, int posNegSwitch, int process);
  double function(double u, double y, int process);

  //Elastic.h//
  RateElastic getRateElasTotal(double p, double T);
  RateElastic getRateElasPos(double u, double T);
  RateElastic getRateElasNeg(double u, double T);

  RateConversion getRateConv(double p, double T);

  double getEnergyTransfer(double p, double T, int process);
  double getMomentumTransfer(double p, double omega, double T, int process);

  double areaOmega(double u, int posNegSwitch, int process);
  double areaQ(double u, double omega, int process);
  double functionOmega(double u, double y, int process) {
    return getElasticRateOmega(u, y, process);}
  double functionQ(double u, double omega, double q, int process) {
    return getElasticRateQ(u, omega, q, process);}
  Jetscape::FourVector getNewMomentumElas(Jetscape::FourVector pVec, double omega, double q);

  //Import.h//
  void readRadiativeRate(Gamma_info *dat, dGammas *Gam);
  void readElasticRateOmega();
  void readElasticRateQ();

  double getRate_qqg(double p, double k);
  double getRate_gqq(double p, double k);
  double getRate_ggg(double p, double k);
  double getRate_qqgamma(double p, double k);
  double use_table(double p, double k, double dGamma[NP][NK], int which_kind);

  double getElasticRateOmega(double u, double omega, int process);
  double getElasticRateQ(double u, double omega, double q, int process);
  double use_elastic_table_omega(double omega, int which_kind);
  double use_elastic_table_q(double u, double omega, int which_kind);

 protected:
  uniform_real_distribution<double> ZeroOneDistribution;

  string PathToTables;
};

double LambertW(double z);

#endif // MARTINI_H

