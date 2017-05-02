#include <fstream>
#include <iomanip>
#include <iomanip>
#include <TF1.h>
#include <TF2.h>
#include <TGraph.h>

#include "fld.h"
#include "eos.h"
#include "icGubser.h"
#include "inc.h"
#include "s95p.h"

using namespace std;

ICGubser::ICGubser() {
}

ICGubser::~ICGubser(void) {}


void ICGubser::setIC(Fluid *f, EoS *eos, double tau) {
  double e, nb, nq, vx = 0., vy = 0., vz = 0.;
  Cell *c;

  double avv_num = 0., avv_den = 0.;
  double Etotal = 0.0;

  for (int ix = 0; ix < f->getNX(); ix++)
    for (int iy = 0; iy < f->getNY(); iy++)
      for (int iz = 0; iz < f->getNZ(); iz++) {
        c = f->getCell(ix, iy, iz);
        double x = f->getX(ix);
        double y = f->getY(iy);
        double eta = f->getZ(iz);

        const double _t = 1.0;
        const double q = 1.0;
        const double r = sqrt(x * x + y * y);
        const double _k =
            2. * q * q * _t * r / (1. + q * q * _t * _t + q * q * r * r);
        e = 4. * q * q / (1. + 2. * q * q * (_t * _t + r * r) +
                          pow(q * q * (_t * _t - r * r), 2)) /
            _t;
        if (e < 0.) e = 0.;
        e = pow(e, 4. / 3.);
        vx = x / (r + 1e-50) * _k;
        vy = y / (r + 1e-50) * _k;
        nb = nq = 0.0;
        vz = 0.0;

        avv_num += sqrt(vx * vx + vy * vy) * e;
        avv_den += e;

        c->setPrimVar(eos, tau, e, nb, nq, 0., vx, vy, vz);
        double _p = eos->p(e, nb, nq, 0.);
        const double gamma2 = 1.0 / (1.0 - vx * vx - vy * vy - vz * vz);
        Etotal +=
            ((e + _p) * gamma2 * (cosh(eta) + vz * sinh(eta)) - _p * cosh(eta));
        c->saveQprev();

        if (e > 0.) c->setAllM(1.);
      }
  cout << "average initial flow = " << avv_num / avv_den << endl;
  cout << "total energy = " << Etotal *f->getDx() * f->getDy() * f->getDz() *
                                   tau << endl;
}
