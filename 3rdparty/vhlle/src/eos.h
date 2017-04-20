#pragma once
#include <cmath>
#include <string>

class TGraph;

// NOTATIONS:
//  e = energy density (local rest frame), [GeV/fm^3]
//  p =  equilibrium pressure, [GeV/fm^3]
// nb =   baryon charge density [1/fm^3]
// nq = electric charge density [1/fm^3]
// ns =  strange charge density [1/fm^3]
// s  = entropy density [1/fm^3]
// T  = temperature  [GeV]
// mub = baryon chemical potential [GeV]
// muq = electric chemical potential [GeV]
// mus = strange chemical potential [GeV]

// abstract EoS class.
// actual EoSes are implemented in derived classes
class EoS {
 public:
  virtual ~EoS() {}
  // eos() gets all EoS relations together:
  // {p,T,mu_b,mu_q,mu_s}={p,T,mu_b,mu_q,mu_s}(e,n_b,n_q,n_s)
  virtual void eos(double e, double nb, double nq, double ns, double &T,
                   double &mub, double &muq, double &mus, double &p) = 0;
  // gets only pressure : p=p(e,n_b,n_q,n_s)
  virtual double p(double e, double nb, double ns, double nq) = 0;
  // gets entropy density
  double s(double e, double nb, double nq, double ns);
  // speed of sound squared: this variant is only used in
  // HLLE solver, where the optimal value is 1/3
  inline double cs2(void) { return 1. / 3.; }
  inline double cs(void) { return sqrt(1. / 3.); }
  // speed of sound squared as a function of energy density
  virtual inline double cs2(double e) { return 1. / 3.; }
  ;
};

// EoS class implementing two variants:
// 1) "SIMPLE": EoS for ultrarelativistic maseless gas
// 2) "TABLE" : EoS p=p(e) from a table
// each variant is enabled by compiling with -D SIMPLE / -D TABLE
class EoSs : public EoS {
 private:
  TGraph *gp, *gT, *gmu;

 public:
  EoSs(std::string fname, int ncols);
  ~EoSs();

  virtual inline void eos(double e, double nb, double nq, double ns, double &T,
                          double &mub, double &muq, double &mus, double &_p) {
    _p = p(e);
    T = t(e);
    mub = muq = mus = 0.;
  }
  virtual inline double p(double e, double nb, double ns, double nq) {
    return p(e);
  }

  double p(double e);
  double dpe(double e);
  double t(double e);
  double mu(double e);

  virtual double cs2(double e) { return dpe(e); }
  // virtual double cs(double e) { return sqrt(dpe(e)) ; }
};
