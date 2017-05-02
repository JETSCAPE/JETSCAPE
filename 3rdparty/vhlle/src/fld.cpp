/******************************************************************************
*                                                                             *
*            vHLLE : a 3D viscous hydrodynamic code                           *
*            version 1.0,            November 2013                            *
*            by Iurii Karpenko                                                *
*  contact:  yu.karpenko@gmail.com                                            *
*  For the detailed description please refer to:                              *
*  http://arxiv.org/abs/1312.4160                                             *
*                                                                             *
*  This code can be freely used and redistributed, provided that this         *
*  copyright appear in all the copies. If you decide to make modifications    *
*  to the code, please contact the authors, especially if you plan to publish *
* the results obtained with such modified code. Any publication of results    *
* obtained using this code must include the reference to                      *
* arXiv:1312.4160 [nucl-th] or the published version of it, when available.   *
*                                                                             *
*******************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "inc.h"
#include "rmn.h"
#include "fld.h"
#include "cll.h"
#include "eos.h"
#include "trancoeff.h"

#define OUTPI

using namespace std;

// returns the velocities in cartesian coordinates, fireball rest frame.
// Y=longitudinal rapidity of fluid
void Fluid::getCMFvariables(Cell *c, double tau, double &e, double &nb,
                            double &nq, double &ns, double &vx, double &vy,
                            double &Y) {
  double p, vz;
  c->getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz);
  double eta = getZ(c->getZ());
  //	Y = eta + TMath::ATanH(vz) ;
  Y = eta + 1. / 2. * log((1. + vz) / (1. - vz));
  vx = vx * cosh(Y - eta) / cosh(Y);
  vy = vy * cosh(Y - eta) / cosh(Y);
}

Fluid::Fluid(EoS *_eos, EoS *_eosH, TransportCoeff *_trcoeff, int _nx, int _ny,
             int _nz, double _minx, double _maxx, double _miny, double _maxy,
             double _minz, double _maxz, double _dt, double eCrit) {
  eos = _eos;
  eosH = _eosH;
  trcoeff = _trcoeff;
  nx = _nx;
  ny = _ny;
  nz = _nz;
  minx = _minx;
  maxx = _maxx;
  miny = _miny;
  maxy = _maxy;
  minz = _minz;
  maxz = _maxz;
  dx = (maxx - minx) / (nx - 1);
  dy = (maxy - miny) / (ny - 1);
  dz = (maxz - minz) / (nz - 1);
  dt = _dt;

  cell = new Cell[nx * ny * nz];

  cell0 = new Cell;
  cell0->setPrimVar(eos, 1.0, 0., 0., 0., 0., 0., 0.,
                    0.);  // tau is not important here, since *0
  cell0->setAllM(0.);

  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++)
      for (int iz = 0; iz < nz; iz++) {
        getCell(ix, iy, iz)->setPrev(X_, getCell(ix - 1, iy, iz));
        getCell(ix, iy, iz)->setNext(X_, getCell(ix + 1, iy, iz));
        getCell(ix, iy, iz)->setPrev(Y_, getCell(ix, iy - 1, iz));
        getCell(ix, iy, iz)->setNext(Y_, getCell(ix, iy + 1, iz));
        getCell(ix, iy, iz)->setPrev(Z_, getCell(ix, iy, iz - 1));
        getCell(ix, iy, iz)->setNext(Z_, getCell(ix, iy, iz + 1));
        getCell(ix, iy, iz)->setPos(ix, iy, iz);
      }

  output_nt = 0;
  output_nx = 0;
  output_ny = 0;
}

Fluid::~Fluid() {
  delete[] cell;
  delete cell0;
}

void Fluid::initOutput(char *dir, int maxstep, double tau0, int cmpr2dOut) {
  //    directory = dir ;
  compress2dOut = cmpr2dOut;
  cout << "maxstep=" << maxstep << endl;
  char command[255];
  sprintf(command, "mkdir -p %s", dir);
  int return_mkdir = system(command);
  cout << "mkdir returns: " << return_mkdir << endl;
  string outx = dir;
  outx.append("/outx.dat");
  string outxvisc = dir;
  outxvisc.append("/outx.visc.dat");
  string outyvisc = dir;
  outyvisc.append("/outy.visc.dat");
  string outdiagvisc = dir;
  outdiagvisc.append("/diag.visc.dat");
  string outy = dir;
  outy.append("/outy.dat");
  string outdiag = dir;
  outdiag.append("/outdiag.dat");
  string outz = dir;
  outz.append("/outz.dat");
  string outaniz = dir;
  outaniz.append("/out.aniz.dat");
  string out2d = dir;
  out2d.append("/out2D.dat");
  string outfreeze = dir;
  outfreeze.append("/freezeout.dat");
  foutx.open(outx.c_str());
  fouty.open(outy.c_str());
  foutz.open(outz.c_str());
  foutdiag.open(outdiag.c_str());
  fout2d.open(out2d.c_str());
  foutxvisc.open(outxvisc.c_str());
  foutyvisc.open(outyvisc.c_str());
  foutdiagvisc.open(outdiagvisc.c_str());
  fout_aniz.open(outaniz.c_str());
  ffreeze.open(outfreeze.c_str());
  //################################################################
  // important remark. for correct diagonal output, nx=ny must hold.
  //################################################################
  foutxvisc << maxstep + 1 << "  " << getNX() << endl;
  foutyvisc << maxstep + 1 << "  " << getNY() << endl;
  foutdiagvisc << maxstep + 1 << "  " << getNX() << endl;
  foutx << "# " << maxstep + 1 << "  " << getNX() << endl;
  fouty << "# " << maxstep + 1 << "  " << getNY() << endl;
  foutdiag << "# " << maxstep + 1 << "  " << getNX() << endl;
  foutz << "# " << maxstep + 1 << "  " << getNZ() << endl;
  fout2d << " " << maxstep + 1 << "  " << (getNX() - 5) + 1 << "  "
         << (getNY() - 5) + 1 << endl;
  fout2d << tau0 << "  " << tau0 + 0.05 * maxstep << "  " << getX(2) << "  "
         << getX(getNX() - 3) << "  " << getY(2) << "  " << getY(getNY() - 3)
         << endl;
  outputGnuplot(tau0);
  fout_aniz << "#  tau  <<v_T>>  e_p  e'_p  (to compare with SongHeinz)\n";
}

void Fluid::correctImagCells(void) {
  double Q[7];
  // Z
  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++) {
      // left boundary
      getCell(ix, iy, 2)->getQ(Q);
      getCell(ix, iy, 1)->setQ(Q);
      getCell(ix, iy, 0)->setQ(Q);
      // right boundary
      getCell(ix, iy, nz - 3)->getQ(Q);
      getCell(ix, iy, nz - 2)->setQ(Q);
      getCell(ix, iy, nz - 1)->setQ(Q);
    }
  // Y
  for (int ix = 0; ix < nx; ix++)
    for (int iz = 0; iz < nz; iz++) {
      // left boundary
      getCell(ix, 2, iz)->getQ(Q);
      getCell(ix, 1, iz)->setQ(Q);
      getCell(ix, 0, iz)->setQ(Q);
      // right boundary
      getCell(ix, ny - 3, iz)->getQ(Q);
      getCell(ix, ny - 2, iz)->setQ(Q);
      getCell(ix, ny - 1, iz)->setQ(Q);
    }
  // X
  for (int iy = 0; iy < ny; iy++)
    for (int iz = 0; iz < nz; iz++) {
      // left boundary
      getCell(2, iy, iz)->getQ(Q);
      getCell(1, iy, iz)->setQ(Q);
      getCell(0, iy, iz)->setQ(Q);
      // right boundary
      getCell(nx - 3, iy, iz)->getQ(Q);
      getCell(nx - 2, iy, iz)->setQ(Q);
      getCell(nx - 1, iy, iz)->setQ(Q);
    }
}

void Fluid::correctImagCellsFull(void) {
  double Q[7], _pi[4][4], _Pi;
  // Z
  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++) {
      // left boundary
      getCell(ix, iy, 2)->getQ(Q);
      getCell(ix, iy, 1)->setQ(Q);
      getCell(ix, iy, 0)->setQ(Q);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) _pi[i][j] = getCell(ix, iy, 2)->getpi(i, j);
      _Pi = getCell(ix, iy, 2)->getPi();

      for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
          getCell(ix, iy, 0)->setpi(i, j, _pi[i][j]);
          getCell(ix, iy, 1)->setpi(i, j, _pi[i][j]);
        }
      getCell(ix, iy, 0)->setPi(_Pi);
      getCell(ix, iy, 1)->setPi(_Pi);
      // right boundary
      getCell(ix, iy, nz - 3)->getQ(Q);
      getCell(ix, iy, nz - 2)->setQ(Q);
      getCell(ix, iy, nz - 1)->setQ(Q);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          _pi[i][j] = getCell(ix, iy, nz - 3)->getpi(i, j);
      _Pi = getCell(ix, iy, nz - 3)->getPi();

      for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
          getCell(ix, iy, nz - 2)->setpi(i, j, _pi[i][j]);
          getCell(ix, iy, nz - 1)->setpi(i, j, _pi[i][j]);
        }
      getCell(ix, iy, nz - 2)->setPi(_Pi);
      getCell(ix, iy, nz - 1)->setPi(_Pi);
    }
  // Y
  for (int ix = 0; ix < nx; ix++)
    for (int iz = 0; iz < nz; iz++) {
      // left boundary
      getCell(ix, 2, iz)->getQ(Q);
      getCell(ix, 1, iz)->setQ(Q);
      getCell(ix, 0, iz)->setQ(Q);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) _pi[i][j] = getCell(ix, 2, iz)->getpi(i, j);
      _Pi = getCell(ix, 2, iz)->getPi();

      for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
          getCell(ix, 0, iz)->setpi(i, j, _pi[i][j]);
          getCell(ix, 1, iz)->setpi(i, j, _pi[i][j]);
        }
      getCell(ix, 0, iz)->setPi(_Pi);
      getCell(ix, 1, iz)->setPi(_Pi);
      // right boundary
      getCell(ix, ny - 3, iz)->getQ(Q);
      getCell(ix, ny - 2, iz)->setQ(Q);
      getCell(ix, ny - 1, iz)->setQ(Q);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          _pi[i][j] = getCell(ix, ny - 3, iz)->getpi(i, j);
      _Pi = getCell(ix, ny - 3, iz)->getPi();

      for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
          getCell(ix, ny - 2, iz)->setpi(i, j, _pi[i][j]);
          getCell(ix, ny - 1, iz)->setpi(i, j, _pi[i][j]);
        }
      getCell(ix, ny - 2, iz)->setPi(_Pi);
      getCell(ix, ny - 1, iz)->setPi(_Pi);
    }
  // X
  for (int iy = 0; iy < ny; iy++)
    for (int iz = 0; iz < nz; iz++) {
      // left boundary
      getCell(2, iy, iz)->getQ(Q);
      getCell(1, iy, iz)->setQ(Q);
      getCell(0, iy, iz)->setQ(Q);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) _pi[i][j] = getCell(2, iy, iz)->getpi(i, j);
      _Pi = getCell(2, iy, iz)->getPi();

      for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
          getCell(0, iy, iz)->setpi(i, j, _pi[i][j]);
          getCell(1, iy, iz)->setpi(i, j, _pi[i][j]);
        }
      getCell(0, iy, iz)->setPi(_Pi);
      getCell(1, iy, iz)->setPi(_Pi);
      // right boundary
      getCell(nx - 3, iy, iz)->getQ(Q);
      getCell(nx - 2, iy, iz)->setQ(Q);
      getCell(nx - 1, iy, iz)->setQ(Q);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          _pi[i][j] = getCell(nx - 3, iy, iz)->getpi(i, j);
      _Pi = getCell(nx - 3, iy, iz)->getPi();

      for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
          getCell(nx - 2, iy, iz)->setpi(i, j, _pi[i][j]);
          getCell(nx - 1, iy, iz)->setpi(i, j, _pi[i][j]);
        }
      getCell(nx - 2, iy, iz)->setPi(_Pi);
      getCell(nx - 1, iy, iz)->setPi(_Pi);
    }
}

void Fluid::updateM(double tau, double dt) {
  for (int ix = 0; ix < getNX(); ix++)
    for (int iy = 0; iy < getNY(); iy++)
      for (int iz = 0; iz < getNZ(); iz++) {
        Cell *c = getCell(ix, iy, iz);
        c->setDM(X_, 0.);
        c->setDM(Y_, 0.);
        c->setDM(Z_, 0.);
        if (getCell(ix, iy, iz)->getMaxM() < 1.) {
          if (getCell(ix + 1, iy, iz)->getM(X_) >= 1. ||
              getCell(ix - 1, iy, iz)->getM(X_) >= 1.)
            c->setDM(X_, dt / dx);
          if (getCell(ix, iy + 1, iz)->getM(Y_) >= 1. ||
              getCell(ix, iy - 1, iz)->getM(Y_) >= 1.)
            c->setDM(Y_, dt / dy);
          if (getCell(ix, iy, iz + 1)->getM(Z_) >= 1. ||
              getCell(ix, iy, iz - 1)->getM(Z_) >= 1.)
            c->setDM(Z_, dt / dz / tau);

          if (c->getDM(X_) == 0. && c->getDM(Y_) == 0.) {
            if (getCell(ix + 1, iy + 1, iz)->getMaxM() >= 1. ||
                getCell(ix + 1, iy - 1, iz)->getMaxM() >= 1. ||
                getCell(ix - 1, iy + 1, iz)->getMaxM() >= 1. ||
                getCell(ix - 1, iy - 1, iz)->getMaxM() >= 1.) {
              c->setDM(X_, 0.707 * dt / dx);
              c->setDM(Y_, 0.707 * dt / dy);
            }
          }
        }  // if
      }

  for (int ix = 0; ix < getNX(); ix++)
    for (int iy = 0; iy < getNY(); iy++)
      for (int iz = 0; iz < getNZ(); iz++) {
        Cell *c = getCell(ix, iy, iz);
        c->addM(X_, c->getDM(X_));
        c->addM(Y_, c->getDM(Y_));
        c->addM(Z_, c->getDM(Z_));
      }
}


void Fluid::outputGnuplot(double tau) {
  double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz;

  // X direction
  for (int ix = 0; ix < nx; ix++) {
    double x = getX(ix);
    Cell *c = getCell(ix, ny / 2, nz / 2);
    getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
    eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
    foutx << setw(14) << tau << setw(14) << x << setw(14) << vx << setw(14)
          << vy << setw(14) << e << setw(14) << nb << setw(14) << t << setw(14)
          << mub;
    foutx << setw(14) << c->getpi(0, 0) << setw(14) << c->getpi(0, 1)
          << setw(14) << c->getpi(0, 2);
    foutx << setw(14) << c->getpi(0, 3) << setw(14) << c->getpi(1, 1)
          << setw(14) << c->getpi(1, 2);
    foutx << setw(14) << c->getpi(1, 3) << setw(14) << c->getpi(2, 2)
          << setw(14) << c->getpi(2, 3);
    foutx << setw(14) << c->getpi(3, 3) << setw(14) << c->getPi() << setw(14)
          << c->getViscCorrCutFlag() << endl;
  }
  foutx << endl;

  // Y direction
  for (int iy = 0; iy < ny; iy++) {
    double y = getY(iy);
    Cell *c = getCell(nx / 2, iy, nz / 2);
    getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
    eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
    fouty << setw(14) << tau << setw(14) << y << setw(14) << vy << setw(14)
          << vx << setw(14) << e << setw(14) << nb << setw(14) << t << setw(14)
          << mub;
    fouty << setw(14) << c->getpi(0, 0) << setw(14) << c->getpi(0, 1)
          << setw(14) << c->getpi(0, 2);
    fouty << setw(14) << c->getpi(0, 3) << setw(14) << c->getpi(1, 1)
          << setw(14) << c->getpi(1, 2);
    fouty << setw(14) << c->getpi(1, 3) << setw(14) << c->getpi(2, 2)
          << setw(14) << c->getpi(2, 3);
    fouty << setw(14) << c->getpi(3, 3) << setw(14) << c->getPi() << setw(14)
          << c->getViscCorrCutFlag() << endl;
  }
  fouty << endl;

  // diagonal
  for (int ix = 0; ix < nx; ix++) {
    double x = getY(ix);
    Cell *c = getCell(ix, ix, nz / 2);
    getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
    eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
    foutdiag << setw(14) << tau << setw(14) << sqrt(2.) * x << setw(14) << vx
             << setw(14) << vy << setw(14) << e << setw(14) << nb << setw(14)
             << t << setw(14) << mub << endl;
    foutdiag << setw(14) << c->getpi(0, 0) << setw(14) << c->getpi(0, 1)
             << setw(14) << c->getpi(0, 2);
    foutdiag << setw(14) << c->getpi(0, 3) << setw(14) << c->getpi(1, 1)
             << setw(14) << c->getpi(1, 2);
    foutdiag << setw(14) << c->getpi(1, 3) << setw(14) << c->getpi(2, 2)
             << setw(14) << c->getpi(2, 3);
    foutdiag << setw(14) << c->getpi(3, 3) << setw(14) << c->getPi() << setw(14)
             << c->getViscCorrCutFlag() << endl;
  }
  foutdiag << endl;

  // Z direction
  for (int iz = 0; iz < nz; iz++) {
    double z = getZ(iz);
    Cell *c = getCell(nx / 2, ny / 2, iz);
    getCMFvariables(getCell(nx / 2, ny / 2, iz), tau, e, nb, nq, ns, vx, vy,
                    vz);
    eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
    foutz << setw(14) << tau << setw(14) << z << setw(14) << vz << setw(14)
          << vx << setw(14) << e << setw(14) << nb << setw(14) << t << setw(14)
          << mub;
    foutz << setw(14) << c->getpi(0, 0) << setw(14) << c->getpi(0, 1)
          << setw(14) << c->getpi(0, 2);
    foutz << setw(14) << c->getpi(0, 3) << setw(14) << c->getpi(1, 1)
          << setw(14) << c->getpi(1, 2);
    foutz << setw(14) << c->getpi(1, 3) << setw(14) << c->getpi(2, 2)
          << setw(14) << c->getpi(2, 3);
    foutz << setw(14) << c->getpi(3, 3) << setw(14) << c->getPi() << setw(14)
          << c->getViscCorrCutFlag() << endl;
  }
  foutz << endl;

  // averaged quantities for y=0
  double T0xx = 0.0, T0yy = 0.0, Txx = 0.0, Tyy = 0.0, vtNum = 0.0, vtDen = 0.0;
  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++) {
      Cell *c = getCell(ix, iy, nz / 2);
      getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
      eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
      double dTxx =
          (e + p) * vx * vx / (1.0 - vx * vx - vy * vy - pow(tanh(vz), 2)) + p;
      double dTyy =
          (e + p) * vy * vy / (1.0 - vx * vx - vy * vy - pow(tanh(vz), 2)) + p;
      T0xx += dTxx;
      T0yy += dTyy;
      Txx += dTxx + c->getpi(1, 1);
      Tyy += dTyy + c->getpi(2, 2);
      vtNum += sqrt(vx * vx + vy * vy) * e /
               sqrt(1.0 - vx * vx - vy * vy - pow(tanh(vz), 2));
      vtDen += e / sqrt(1.0 - vx * vx - vy * vy - pow(tanh(vz), 2));
    }
  fout_aniz << setw(14) << tau << setw(14) << vtNum / vtDen << setw(14)
            << (T0xx - T0yy) / (T0xx + T0yy) << setw(14)
            << (Txx - Tyy) / (Txx + Tyy) << endl;
}

// unput: geom. rapidity + velocities in Bjorken frame, --> output: velocities
// in lab.frame
void transformToLab(double eta, double &vx, double &vy, double &vz) {
  const double Y = eta + 1. / 2. * log((1. + vz) / (1. - vz));
  vx = vx * cosh(Y - eta) / cosh(Y);
  vy = vy * cosh(Y - eta) / cosh(Y);
  vz = tanh(Y);
}

void Fluid::calcTotals(double tau) {
  double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, Q[7];
  double E = 0., Efull = 0., Px = 0., Nb1 = 0., Nb2 = 0., S = 0.;
  double eta = 0;

  fout2d << endl;
  for (int ix = 2; ix < nx - 2; ix++)
    for (int iy = 2; iy < ny - 2; iy++)
      for (int iz = 2; iz < nz - 2; iz++) {
        Cell *c = getCell(ix, iy, iz);
        getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
        c->getQ(Q);
        eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
        double s = eos->s(e, nb, nq, ns);
        eta = getZ(iz);
        E += tau * (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                 (cosh(eta) - tanh(vz) * sinh(eta)) -
             tau * p * cosh(eta);
        Nb1 += Q[NB_];
        Nb2 += tau * nb * (cosh(eta) - tanh(vz) * sinh(eta)) /
               sqrt(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
        // -- noneq. corrections to entropy flux
        const double gmumu[4] = {1., -1., -1., -1.};
        double deltas = 0.;
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++)
            deltas += pow(c->getpi(i, j), 2) * gmumu[i] * gmumu[j];
        if (t > 0.05) s += 1.5 * deltas / ((e + p) * t);
        S += tau * s * (cosh(eta) - tanh(vz) * sinh(eta)) /
             sqrt(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
        Efull +=
            tau * (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                (cosh(eta) - tanh(vz) * sinh(eta)) -
            tau * p * cosh(eta);
        if (trcoeff->isViscous())
          Efull += tau * c->getpi(0, 0) * cosh(eta) +
                   tau * c->getpi(0, 3) * sinh(eta);
      }
  E = E * dx * dy * dz;
  Efull = Efull * dx * dy * dz;
  Nb1 *= dx * dy * dz;
  Nb2 *= dx * dy * dz;
  S *= dx * dy * dz;
  cout << setw(16) << "calcTotals: E = " << setw(14) << E
       << "  Efull = " << setw(14) << Efull << endl;
  cout << setw(16) << "Px = " << setw(14) << Px << "      S = " << setw(14) << S
       << endl;
}
