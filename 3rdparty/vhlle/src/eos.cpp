#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "eos.h"
#include "inc.h"
#include <TGraph.h>

using namespace std;

const double bagp = pow(247.19 / 197.32, 4) / gevtofm;
const double bagt = pow(247.19 / 197.32, 4) / gevtofm;
const double deg = 16.0 + 3.0 * 12.0 * (7.0 / 8.0);

// EoS choise --> MOVED to Makefile
//#define TABLE // Laine, etc
//#define SIMPLE  // p=e/3

double EoS::s(double e, double nb, double nq, double ns) {
  double T, mub, muq, mus, p;
  eos(e, nb, nq, ns, T, mub, muq, mus, p);
  if (T > 0.0)
    return (e + p - mub * nb - muq * nq - mus * ns) / T;
  else
    return 0.;
}

EoSs::EoSs(string fname, int ncols) {
#if defined TABLE || defined LAINE_CFO

  int edat = 10000;
  double* e = new double[edat];
  double* pGrid = new double[edat];
  double* tpGrid = new double[edat];
  double* muGrid = new double[edat];

  ifstream finput(fname.c_str(), ios::in);
  // ofstream fout("debug.txt");
  if (!finput) {
    cerr << "can't open input file \"" << fname.c_str() << "\"" << endl;
    exit(1);
  }
  edat = 0;
  while (!finput.eof()) {
    if (ncols == 3) {
      finput >> e[edat] >> pGrid[edat] >> tpGrid[edat];
      muGrid[edat] = 0.;
    } else {
      finput >> e[edat] >> pGrid[edat] >> tpGrid[edat] >> muGrid[edat];
    }
    if (pGrid[edat] < 0.) pGrid[edat] = 0.;
    edat++;
  }
  finput.close();

  gp = new TGraph(edat, e, pGrid);
  gT = new TGraph(edat, e, tpGrid);
  gmu = new TGraph(edat, e, muGrid);

#elif defined SIMPLE
// nothing
#endif
}

EoSs::~EoSs(void) {}

double EoSs::p(double e) {
#if defined TABLE
  return gp->Eval(e);
#elif defined SIMPLE
  return e / 3.;
#endif
}

double EoSs::dpe(double e) {
#if defined TABLE
  return (gp->Eval(e * 1.1) - gp->Eval(e)) / (0.1 * e);
#elif defined SIMPLE
  return 1. / 3.;
#endif
}

double EoSs::t(double e) {
#if defined TABLE
  return gT->Eval(e);
#elif defined SIMPLE
  const double cnst =
      (16 + 0.5 * 21.0 * 2.5) * pow(C_PI, 2) / 30.0 / pow(0.197326968, 3);
  return e > 0. ? 1.0 * pow(e / cnst, 0.25) : 0.;
#endif
}

double EoSs::mu(double e) {
#if defined TABLE
  return gmu->Eval(e);
#elif defined SIMPLE
  return 0.;
#endif
}
