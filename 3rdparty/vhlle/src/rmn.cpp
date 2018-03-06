#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <execinfo.h>
#include <signal.h>
#include "rmn.h"

using namespace std;

// flag optionally needed to track
// problems in the transformation procedures
bool debugRiemann;

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  cout << "Error: signal " << sig << endl;
  backtrace_symbols_fd(array, size, 2);
  exit(1);
}

void transformPV(EoS *eos, double Q[7], double &e, double &p, double &nb,
                 double &nq, double &ns, double &vx, double &vy, double &vz) {
  // conserved -> primitive transtormation requires
  // a numerical solution to 1D nonlinear algebraic equation:
  // v = M / ( Q_t + p(Q_t-M*v, n) )       (A.2)
  // M being the modulo of the vector {Q_x, Q_y, Q_z}.
  // Bisection/Newton methods are used to solve the equation.
  const int MAXIT = 100;        // maximum number of iterations
  const double dpe = 1. / 3.;   // dp/de estimate for Newton method
  const double corrf = 0.9999;  // corrected value of M
  // when it brakes the speed of light limit, M>Q_t
  double v, vl = 0., vh = 1., dvold, dv, f, df;
  if (debugRiemann) {
    cout << "transformPV debug---------------\n";
    cout << setw(14) << Q[0] << setw(14) << Q[1] << setw(14) << Q[2] << setw(14)
         << Q[3] << endl;
    cout << setw(14) << Q[4] << setw(14) << Q[5] << setw(14) << Q[6] << endl;
  }
  double M = sqrt(Q[X_] * Q[X_] + Q[Y_] * Q[Y_] + Q[Z_] * Q[Z_]);
  if (Q[T_] <= 0.) {
    e = 0.;
    p = 0.;
    vx = vy = vz = 0.;
    nb = nq = ns = 0.;
    return;
  }
  if (M == 0.) {
    e = Q[T_];
    vx = 0.;
    vy = 0.;
    vz = 0.;
    nb = Q[NB_];
    nq = Q[NQ_];
    ns = Q[NS_];
    p = eos->p(e, nb, nq, ns);
    return;
  }
  if (M > Q[T_]) {
    Q[X_] *= corrf * Q[T_] / M;
    Q[Y_] *= corrf * Q[T_] / M;
    Q[Z_] *= corrf * Q[T_] / M;
    M = Q[T_] * corrf;
  }

  v = 0.5 * (vl + vh);
  e = Q[T_] - M * v;
  if (e < 0.) e = 0.;
  nb = Q[NB_] * sqrt(1 - v * v);
  nq = Q[NQ_] * sqrt(1 - v * v);
  ns = Q[NS_] * sqrt(1 - v * v);
  p = eos->p(e, nb, nq, ns);
  f = (Q[T_] + p) * v - M;
  df = (Q[T_] + p) - M * v * dpe;
  dvold = vh - vl;
  dv = dvold;
  for (int i = 0; i < MAXIT; i++) {
    if ((f + df * (vh - v)) * (f + df * (vl - v)) > 0. ||
        fabs(2. * f) > fabs(dvold * df)) {  // bisection
      dvold = dv;
      dv = 0.5 * (vh - vl);
      v = vl + dv;
      //			cout << "BISECTION v = " << setw(12) << v << endl ;
    } else {  // Newton
      dvold = dv;
      dv = f / df;
      v -= dv;
      //			cout << "NEWTON v = " << setw(12) << v << endl ;
    }
    if (fabs(dv) < 0.00001) break;

    e = Q[T_] - M * v;
    if (e < 0.) e = 0.;
    nb = Q[NB_] * sqrt(1 - v * v);
    nq = Q[NQ_] * sqrt(1 - v * v);
    ns = Q[NS_] * sqrt(1 - v * v);
    p = eos->p(e, nb, nq, ns);
    f = (Q[T_] + p) * v - M;
    df = (Q[T_] + p) - M * v * dpe;

    if (f > 0.)
      vh = v;
    else
      vl = v;
    if (debugRiemann)
      cout << "step " << i << "  " << e << "  " << nb << "  " << nq << "  "
           << ns << "  " << p << "  " << vx << "  " << vy << "  " << vz << endl;
    //		if(i==40) { cout << "error : does not converge\n" ; exit(1) ; } ;
  }  // for loop
     //----------after
     //	v = 0.5*(vh+vl) ;
  vx = v * Q[X_] / M;
  vy = v * Q[Y_] / M;
  vz = v * Q[Z_] / M;
  e = Q[T_] - M * v;
  p = eos->p(e, nb, nq, ns);
  nb = Q[NB_] * sqrt(1 - vx * vx - vy * vy - vz * vz);
  nq = Q[NQ_] * sqrt(1 - vx * vx - vy * vy - vz * vz);
  ns = Q[NS_] * sqrt(1 - vx * vx - vy * vy - vz * vz);
  if (e < 0. || sqrt(vx * vx + vy * vy + vz * vz) > 1.) {
    cout << Q[T_] << "  " << Q[X_] << "  " << Q[Y_] << "  " << Q[Z_] << "  "
         << Q[NB_] << endl;
    cout << "transformRF::Error\n";
  }
  if (!(nb < 0. || nb >= 0.)) {
    cout << "transformRF::Error nb=#ind\n";
    //		return ;
  }
}

void transformPVBulk(EoS *eos, double Pi, double Q[7], double &e, double &p,
                     double &nb, double &nq, double &ns, double &vx, double &vy,
                     double &vz) {
  const int MAXIT = 100;
  const double dpe = 1. / 3.;
  const double corrf = 0.9999;
  double v, vl = 0., vh = 1., dvold, dv, f, df;
  if (debugRiemann) {
    cout << "transformPVBulk debug---------------\n";
    cout << setw(14) << Q[0] << setw(14) << Q[1] << setw(14) << Q[2] << setw(14)
         << Q[3] << endl;
    cout << setw(14) << Q[4] << setw(14) << Q[5] << setw(14) << Q[6] << endl;
  }
  double M = sqrt(Q[X_] * Q[X_] + Q[Y_] * Q[Y_] + Q[Z_] * Q[Z_]);
  if (Q[T_] <= 0.) {
    e = 0.;
    p = 0.;
    vx = vy = vz = 0.;
    nb = nq = ns = 0.;
    return;
  }
  if (M == 0.) {
    e = Q[T_];
    vx = 0.;
    vy = 0.;
    vz = 0.;
    nb = Q[NB_];
    nq = Q[NQ_];
    ns = Q[NS_];
    p = eos->p(e, nb, nq, ns) + Pi;
    return;
  }
  if (M > Q[T_]) {
    Q[X_] *= corrf * Q[T_] / M;
    Q[Y_] *= corrf * Q[T_] / M;
    Q[Z_] *= corrf * Q[T_] / M;
    M = Q[T_] * corrf;
  }

  v = 0.5 * (vl + vh);
  e = Q[T_] - M * v;
  if (e < 0.) e = 0.;
  nb = Q[NB_] * sqrt(1 - v * v);
  nq = Q[NQ_] * sqrt(1 - v * v);
  ns = Q[NS_] * sqrt(1 - v * v);
  p = eos->p(e, nb, nq, ns) + Pi;
  f = (Q[T_] + p) * v - M;
  df = (Q[T_] + p) - M * v * dpe;
  dvold = vh - vl;
  dv = dvold;
  for (int i = 0; i < MAXIT; i++) {
    if ((f + df * (vh - v)) * (f + df * (vl - v)) > 0. ||
        fabs(2. * f) > fabs(dvold * df)) {  // bisection
      dvold = dv;
      dv = 0.5 * (vh - vl);
      v = vl + dv;
      //			cout << "BISECTION v = " << setw(12) << v << endl ;
    } else {  // Newton
      dvold = dv;
      dv = f / df;
      v -= dv;
      //			cout << "NEWTON v = " << setw(12) << v << endl ;
    }
    if (fabs(dv) < 0.00001) break;

    e = Q[T_] - M * v;
    if (e < 0.) e = 0.;
    nb = Q[NB_] * sqrt(1 - v * v);
    nq = Q[NQ_] * sqrt(1 - v * v);
    ns = Q[NS_] * sqrt(1 - v * v);
    p = eos->p(e, nb, nq, ns) + Pi;
    f = (Q[T_] + p) * v - M;
    df = (Q[T_] + p) - M * v * dpe;

    if (f > 0.)
      vh = v;
    else
      vl = v;
    if (debugRiemann)
      cout << "step " << i << "  " << e << "  " << nb << "  " << nq << "  "
           << ns << "  " << p << "  " << vx << "  " << vy << "  " << vz << endl;
    //		if(i==40) { cout << "error : does not converge\n" ; exit(1) ; } ;
  }  // for loop
     //----------after
     //	v = 0.5*(vh+vl) ;
  vx = v * Q[X_] / M;
  vy = v * Q[Y_] / M;
  vz = v * Q[Z_] / M;
  e = Q[T_] - M * v;
  p = eos->p(e, nb, nq, ns);
  nb = Q[NB_] * sqrt(1 - vx * vx - vy * vy - vz * vz);
  nq = Q[NQ_] * sqrt(1 - vx * vx - vy * vy - vz * vz);
  ns = Q[NS_] * sqrt(1 - vx * vx - vy * vy - vz * vz);
  if (e < 0. || sqrt(vx * vx + vy * vy + vz * vz) > 1.) {
    cout << Q[T_] << "  " << Q[X_] << "  " << Q[Y_] << "  " << Q[Z_] << "  "
         << Q[NB_] << endl;
    cout << "transformRF::Error\n";
  }
  if (!(nb < 0. || nb >= 0.)) {
    cout << "transformRF::Error nb=#ind\n";
    cout << "Q [7]: " << Q[0] << " " << Q[1] << " " << Q[2] << " " << Q[3]
         << " " << Q[4] << " " << Q[5] << " " << Q[6] << endl;
    cout << "e vx vy vz nb nq ns: " << e << " " << vx << " " << vy << " " << vz
         << " " << nb << " " << nq << " " << ns << endl;
    handler(333);
    //		return ;
  }
}

void transformCV(double e, double p, double nb, double nq, double ns, double vx,
                 double vy, double vz, double Q[]) {
  const double gamma2 = 1. / (1 - vx * vx - vy * vy - vz * vz);
  Q[T_] = (e + p * (vx * vx + vy * vy + vz * vz)) * gamma2;
  Q[X_] = (e + p) * vx * gamma2;
  Q[Y_] = (e + p) * vy * gamma2;
  Q[Z_] = (e + p) * vz * gamma2;
  Q[NB_] = nb * sqrt(gamma2);
  Q[NQ_] = nq * sqrt(gamma2);
  Q[NS_] = ns * sqrt(gamma2);
}
