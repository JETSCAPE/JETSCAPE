#include <fstream>
#include <iomanip>
#include <iomanip>
#include <TF1.h>
#include <TF2.h>
#include <TGraph.h>

#include "fld.h"
#include "eos.h"
#include "ic.h"
#include "inc.h"
#include "s95p.h"

using namespace std;


IC::IC(char* icInputFile, double s0ScaleFactor)
{
 s95p::loadSongIC(icInputFile, s0ScaleFactor);
}

IC::~IC(void) {}


void IC::setIC(Fluid *f, EoS *eos, double tau) {
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

        double eta1 = fabs(eta) < 1.3 ? 0.0 : fabs(eta) - 1.3;
        double etaFactor =
              exp(-eta1 * eta1 / 2.1 / 2.1) * (fabs(eta) < 5.3 ? 1.0 : 0.0);
        e = s95p::getSongEps(x, y) * etaFactor;
        if (e < 0.5) e = 0.0;
        vx = vy = 0.0;
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
