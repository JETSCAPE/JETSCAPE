#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

extern int glauberVariable;

namespace s95p {
using namespace std;

int NX, NY;  // N*N table
double** e, xmin, xmax, ymin, ymax;

double s95p_s(double e) {
  double val;
  if (e < 0.1270769021427449)
    val = 12.2304 * pow(e, 1.16849);
  else if (e < 0.4467079524674040)
    val = 11.9279 * pow(e, 1.15635);
  else if (e < 1.9402832534193788)
    val = 0.0580578 + 11.833 * pow(e, 1.16187);
  else if (e < 3.7292474570977285)
    val =
        18.202 * e - 62.021814 - 4.85479 * exp(-2.72407e-11 * pow(e, 4.54886)) +
        65.1272 * pow(e, -0.128012) * exp(-0.00369624 * pow(e, 1.18735)) -
        4.75253 * pow(e, -1.18423);
  else
    val = 18.202 * e - 63.0218 - 4.85479 * exp(-2.72407e-11 * pow(e, 4.54886)) +
          65.1272 * pow(e, -0.128012) * exp(-0.00369624 * pow(e, 1.18735));
  return pow(val, 3. / 4.);
}

double s95p_e(double s) {
  double e = pow(s, 4. / 3.) / 18.2, e0;
  do {
    e0 = e;
    e = e0 - (s95p_s(e0) - s) / (6.608681233 * pow(e0, -0.25));
  } while (fabs(e - e0) > 1e-10);
  return e;
}

void loadSongIC(char* filename, double factor) {
  ifstream fin(filename);
  if (!fin) {
    cout << "cannot open " << filename << endl;
    exit(1);
  }
  fin >> NX >> NY;
  double* x = new double[NX];
  double* y = new double[NY];
  e = new double* [NX];
  for (int ix = 0; ix < NX; ix++) e[ix] = new double[NY];

  for (int ix = 0; ix < NX; ix++)
    for (int iy = 0; iy < NY; iy++) {
      fin >> x[ix] >> y[iy] >> e[ix][iy];
      if (glauberVariable == 1)
        e[ix][iy] = s95p_e(factor * e[ix][iy]);
      else
        e[ix][iy] = factor * e[ix][iy];
    }
  xmin = x[0];
  xmax = x[NX - 1];
  ymin = y[0];
  ymax = y[NY - 1];
  delete[] x;
  delete[] y;
}

double getSongEps(double x, double y) {
  // linear interpolation in 2D is used
  const double dx = (xmax - xmin) / (NX - 1);
  const double dy = (ymax - ymin) / (NY - 1);
  const int ix = (int)((x - xmin) / dx);
  const int iy = (int)((y - ymin) / dy);
  if (ix < 0 || ix > NX - 2 || iy < 0 || iy > NY - 2) return 0.0;
  const double sx = x - xmin - ix * dx;
  const double sy = y - ymin - iy * dy;
  double wx[2] = {(1. - sx / dx), sx / dx};
  double wy[2] = {(1. - sy / dy), sy / dy};
  double result = 0;
  for (int jx = 0; jx < 2; jx++)
    for (int jy = 0; jy < 2; jy++) {
      result += wx[jx] * wy[jy] * e[ix + jx][iy + jy];
    }
  return result;
}

}  // end s95p
