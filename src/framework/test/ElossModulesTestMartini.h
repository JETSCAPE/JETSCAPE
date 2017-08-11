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

#include <fstream>
#include <math.h>
#include <gsl/gsl_sf_lambert.h>
#include "JetEnergyLossModule.h"
#include "constants.h"

//Random.h//
#define NNM 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL         /* Least significant 31 bits */

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
  
 private:

};

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

  double *dGamma_qq;
  double *dGamma_qg;
  double *dGamma_qq_q;
  double *dGamma_qg_q;

  //Random.h//
  /* The array for the state vector */
  unsigned long long mt[NNM]; 
  /* mti==NNM+1 means mt[NNM] is not initialized */
  int mti;

 public:
  
  Martini();
  virtual ~Martini();

  //main//
  void Init();
  void DoEnergyLoss(double deltaT, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
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
  FourVector getNewMomentumElas(FourVector pVec, double omega, double q);

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

  //Random.h//
  void init_genrand64(unsigned long long seed);
  void init_by_array64(unsigned long long init_key[],
        	       unsigned long long key_length);
  unsigned long long genrand64_int64(void);
  long long genrand64_int63(void);
  double genrand64_real1(void);
  double genrand64_real2(void);
  double genrand64_real3(void);

};

#endif

