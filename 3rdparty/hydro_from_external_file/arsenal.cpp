// Ver 1.7.0
// Zhi Qiu
/*==========================================================================================
Change logs: see arsenal.h
==========================================================================================*/

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <cstdarg>

#include "arsenal.h"

#define OUTPUT_PRECISION 10

// for log_gamma function
#define D1      -0.5772156649015328605195174
#define D2      0.4227843350984671393993777
#define D4      1.791759469228055000094023
#define SQRTPI  0.9189385332046727417803297
#define FRTBIG  1.42E+09
#define PNT68   0.6796875
#define XBIG    4.08E+36
#define MACHINE_EPSILON 2.22044604925031e-016

using namespace std;

//**********************************************************************
double sixPoint2dInterp(double x, double y,
    double v00, double v01, double v02, double v10, double v11, double v20)
{
  /* Assume a 2d function f(x,y) has a quadratic form:
    f(x,y) = axx*x^2 + axy*x*y + ayy*y^2 + bx*x + by*y + c
    Also assume that:
    f(0,0)=v00, f(0,1)=v01, f(0,2)=v02, f(1,0)=v10, f(1,1)=v11, f(2,0)=v20
    Then its value at (x,y) can be solved.
    The best result is produced if x and y are between 0 and 1.
  */

  // Calculate coefficients:
  double axx = 1.0/2.0*(v00 - 2*v10 + v20);
  double axy = v00 - v01 - v10 + v11;
  double ayy = 1.0/2.0*(v00 - 2*v01 + v02);
  double bx = 1.0/2.0*(-3.0*v00 + 4*v10 - v20);
  double by = 1.0/2.0*(-3.0*v00 + 4*v01 - v02);
  double c = v00;

  // Calcualte f(x,y):
  return axx*x*x + axy*x*y + ayy*y*y + bx*x + by*y + c;
}


//**********************************************************************
double interpCubicDirect(vector<double>* x, vector<double>* y, double x0)
// Returns the interpreted value of y=y(x) at x=x0 using cubic polynomial interpolation method.
// -- x,y: the independent and dependent double x0ables; x is assumed to be equal spaced and increasing
// -- x0: where the interpolation should be performed
{
  long size = x->size();
  if (size==1) {cout<<"interpCubicDirect warning: table size = 1"; return (*y)[0];}
  double dx = (*x)[1]-(*x)[0]; // increment in x

  // if close to left end:
  if (abs(x0-(*x)[0])<dx*1e-30) return (*y)[0];

  // find x's integer index
  long idx = floor((x0-(*x)[0])/dx);

  if (idx<0 || idx>=size-1)
  {
    cout    << "interpCubicDirect: x0 out of bounds." << endl
            << "x ranges from " << (*x)[0] << " to " << (*x)[size-1] << ", "
            << "x0=" << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl;
    exit(-1);
  }

  if (idx==0)
  {
    // use quadratic interpolation at left end
    double A0 = (*y)[0], A1 = (*y)[1], A2 = (*y)[2], deltaX = x0 - (*x)[0]; // deltaX is the increment of x0 compared to the closest lattice point
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else if (idx==size-2)
  {
    // use quadratic interpolation at right end
    double A0 = (*y)[size-3], A1 = (*y)[size-2], A2 = (*y)[size-1], deltaX = x0 - ((*x)[0] + (idx-1)*dx);
    return (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX - (3.0*A0-4.0*A1+A2)/(2.0*dx)*deltaX + A0;
  }
  else
  {
    // use cubic interpolation
    double A0 = (*y)[idx-1], A1 = (*y)[idx], A2 = (*y)[idx+1], A3 = (*y)[idx+2], deltaX = x0 - ((*x)[0] + idx*dx);
    //cout << A0 << "  " << A1 << "  " << A2 << "  " << A3 << endl;
    return (-A0+3.0*A1-3.0*A2+A3)/(6.0*dx*dx*dx)*deltaX*deltaX*deltaX
            + (A0-2.0*A1+A2)/(2.0*dx*dx)*deltaX*deltaX
            - (2.0*A0+3.0*A1-6.0*A2+A3)/(6.0*dx)*deltaX
            + A1;
  }

}




//**********************************************************************
double interpLinearDirect(vector<double>* x, vector<double>* y, double x0)
// Returns the interpreted value of y=y(x) at x=x0 using linear interpolation method.
// -- x,y: the independent and dependent double x0ables; x is assumed to be equal spaced and increasing
// -- x0: where the interpolation should be performed
{
  long size = x->size();
  if (size==1) {cout<<"interpLinearDirect warning: table size = 1"<<endl; return (*y)[0];}
  double dx = (*x)[1]-(*x)[0]; // increment in x

  // if close to left end:
  if (abs(x0-(*x)[0])<dx*1e-30) return (*y)[0];

  // find x's integer index
  long idx = floor((x0-(*x)[0])/dx);

  if (idx<0 || idx>=size-1)
  {
    cout    << "interpLinearDirect: x0 out of bounds." << endl
            << "x ranges from " << (*x)[0] << " to " << (*x)[size-1] << ", "
            << "x0=" << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl;
    exit(-1);
  }

  return (*y)[idx] + ((*y)[idx+1]-(*y)[idx])/dx*(x0-(*x)[idx]);

}




//**********************************************************************
double interpNearestDirect(vector<double>* x, vector<double>* y, double x0)
// Returns the interpreted value of y=y(x) at x=x0 using nearest interpolation method.
// -- x,y: the independent and dependent double x0ables; x is assumed to be equal spaced and increasing
// -- x0: where the interpolation should be performed
{
  long size = x->size();
  if (size==1) {cout<<"interpNearestDirect warning: table size = 1"<<endl; return (*y)[0];}
  double dx = (*x)[1]-(*x)[0]; // increment in x

  // if close to left end:
  if (abs(x0-(*x)[0])<dx*1e-30) return (*y)[0];

  // find x's integer index
  long idx = floor((x0-(*x)[0])/dx);

  if (idx<0 || idx>=size-1)
  {
    cout    << "interpNearestDirect: x0 out of bounds." << endl
            << "x ranges from " << (*x)[0] << " to " << (*x)[size-1] << ", "
            << "x0=" << x0 << ", " << "dx=" << dx << ", " << "idx=" << idx << endl;
    exit(-1);
  }

  return x0-(*x)[idx]>dx/2 ? (*y)[idx+1] : (*y)[idx];

}




//**********************************************************************
double interpCubicMono(vector<double>* x, vector<double>* y, double xx)
// Note that this function does NOT perform well with small x and y table spacing; in which case use "direct" version instead.
// Returns the interpreted value of y=y(x) at x=x0 using cubic polynomial interpolation method.
// -- x,y: the independent and dependent double x0ables; x is *NOT* assumed to be equal spaced but it has to be increasing
// -- xx: where the interpolation should be performed
{
  long size = x->size();
  if (size==1) {cout<<"interpCubicMono warning: table size = 1"<<endl; return (*y)[0];}

  // if close to left end:
  if (abs(xx-(*x)[0])<((*x)[1]-(*x)[0])*1e-30) return (*y)[0];

  // find x's integer index
  long idx = binarySearch(x, xx);

  if (idx<0 || idx>=size-1)
  {
    cout    << "interpCubicMono: x0 out of bounds." << endl
            << "x ranges from " << (*x)[0] << " to " << (*x)[size-1] << ", "
            << "xx=" << xx << ", " << "idx=" << idx << endl;
    exit(-1);
  }

  if (idx==0)
  {
    // use linear interpolation at the left end
    return (*y)[0] + ( (*y)[1]-(*y)[0] )/( (*x)[1]-(*x)[0] )*( xx-(*x)[0] );
  }
  else if (idx==size-2)
  {
    // use linear interpolation at the right end
    return (*y)[size-2] + ( (*y)[size-1]-(*y)[size-2] )/( (*x)[size-1]-(*x)[size-2] )*( xx-(*x)[size-2] );
  }
  else
  {
    // use cubic interpolation
    long double y0 = (*y)[idx-1], y1 = (*y)[idx], y2 = (*y)[idx+1], y3 = (*y)[idx+2];
    long double y01=y0-y1, y02=y0-y2, y03=y0-y3, y12=y1-y2, y13=y1-y3, y23=y2-y3;
    long double x0 = (*x)[idx-1], x1 = (*x)[idx], x2 = (*x)[idx+1], x3 = (*x)[idx+2];
    long double x01=x0-x1, x02=x0-x2, x03=x0-x3, x12=x1-x2, x13=x1-x3, x23=x2-x3;
    long double x0s=x0*x0, x1s=x1*x1, x2s=x2*x2, x3s=x3*x3;
    long double denominator = x01*x02*x12*x03*x13*x23;
    long double C0, C1, C2, C3;
    C0 = (x0*x02*x2*x03*x23*x3*y1
          + x1*x1s*(x0*x03*x3*y2+x2s*(-x3*y0+x0*y3)+x2*(x3s*y0-x0s*y3))
          + x1*(x0s*x03*x3s*y2+x2*x2s*(-x3s*y0+x0s*y3)+x2s*(x3*x3s*y0-x0*x0s*y3))
          + x1s*(x0*x3*(-x0s+x3s)*y2+x2*x2s*(x3*y0-x0*y3)+x2*(-x3*x3s*y0+x0*x0s*y3))
          )/denominator;
    C1 = (x0s*x03*x3s*y12
          + x2*x2s*(x3s*y01+x0s*y13)
          + x1s*(x3*x3s*y02+x0*x0s*y23-x2*x2s*y03)
          + x2s*(-x3*x3s*y01-x0*x0s*y13)
          + x1*x1s*(-x3s*y02+x2s*y03-x0s*y23)
          )/denominator;
    C2 = (-x0*x3*(x0s-x3s)*y12
          + x2*(x3*x3s*y01+x0*x0s*y13)
          + x1*x1s*(x3*y02+x0*y23-x2*y03)
          + x2*x2s*(-x3*y01-x0*y13)
          + x1*(-x3*x3s*y02+x2*x2s*y03-x0*x0s*y23)
          )/denominator;
    C3 = (x0*x03*x3*y12
          + x2s*(x3*y01+x0*y13)
          + x1*(x3s*y02+x0s*y23-x2s*y03)
          + x2*(-x3s*y01-x0s*y13)
          + x1s*(-x3*y02+x2*y03-x0*y23)
          )/denominator;
/*    cout  << x0s*x03*x3s*y12 << "  "
          <<  x2*x2s*(x3s*y01+x0s*y13) << "   "
          <<  x1s*(x3*x3s*y02+x0*x0s*y23-x2*x2s*y03) << "  "
          <<  x2s*(-x3*x3s*y01-x0*x0s*y13) << "  "
          <<  x1*x1s*(-x3s*y02+x2s*y03-x0s*y23) << endl;
    cout << denominator << endl;

    cout << x0 << " " << x1 << "  " << x2 << "  " << x3 << endl;
    cout << y0 << " " << y1 << "  " << y2 << "  " << y3 << endl;
    cout << C0 << "  " << C1 << "  " << C2 << "  " << C3 << endl;*/
    return C0 + C1*xx + C2*xx*xx + C3*xx*xx*xx;
  }

}





//**********************************************************************
double interpLinearMono(vector<double>* x, vector<double>* y, double xx)
// Returns the interpreted value of y=y(x) at x=x0 using linear interpolation method.
// -- x,y: the independent and dependent double x0ables; x is *NOT* assumed to be equal spaced but it has to be increasing
// -- xx: where the interpolation should be performed
{
  long size = x->size();
  if (size==1) {cout<<"interpLinearMono warning: table size = 1"<<endl; return (*y)[0];}

  // if close to left end:
  if (abs(xx-(*x)[0])<((*x)[1]-(*x)[0])*1e-30) return (*y)[0];

  // find x's integer index
  long idx = binarySearch(x, xx);

  if (idx<0 || idx>=size-1)
  {
    cout    << "interpLinearMono: x0 out of bounds." << endl
            << "x ranges from " << (*x)[0] << " to " << (*x)[size-1] << ", "
            << "xx=" << xx << ", " << "idx=" << idx << endl;
    exit(-1);
  }

  return (*y)[idx] + ( (*y)[idx+1]-(*y)[idx] )/( (*x)[idx+1]-(*x)[idx] )*( xx-(*x)[idx] );

}




//**********************************************************************
double interpNearestMono(vector<double>* x, vector<double>* y, double xx)
// Returns the interpreted value of y=y(x) at x=x0 using nearest interpolation method.
// -- x,y: the independent and dependent double x0ables; x is *NOT* assumed to be equal spaced but it has to be increasing
// -- xx: where the interpolation should be performed
{
  long size = x->size();
  if (size==1) {cout<<"interpNearestMono warning: table size = 1"<<endl; return (*y)[0];}

  // if close to left end:
  if (abs(xx-(*x)[0])<((*x)[1]-(*x)[0])*1e-30) return (*y)[0];

  // find x's integer index
  long idx = binarySearch(x, xx);

  if (idx<0 || idx>=size-1)
  {
    cout    << "interpNearestMono: x0 out of bounds." << endl
            << "x ranges from " << (*x)[0] << " to " << (*x)[size-1] << ", "
            << "xx=" << xx << ", " << "idx=" << idx << endl;
    exit(-1);
  }

  return xx-(*x)[idx] > (*x)[idx+1]-xx ? (*y)[idx+1] : (*y)[idx];

}




//**********************************************************************
double invertFunc(double (*func)(double), double y, double xL, double xR, double dx, double x0, double relative_accuracy)
//Purpose:
//  Return x=func^(-1)(y) using Newton method.
//  -- func: double 1-argument function to be inverted
//  -- xL: left boundary (for numeric derivative)
//  -- xR: right boundary (for numeric derivative)
//  -- dx: step (for numeric derivative)
//  -- x0: initial value
//  -- y: the value to be inverted
//  -- Returns inverted value
//Solve: f(x)=0 with f(x)=table(x)-y => f'(x)=table'(x)
{
  double accuracy;
  int tolerance;

  double XX1, XX2; // used in iterations
  double F0, F1, F2, F3, X1, X2; // intermedia variables
  int impatience; // number of iterations


  // initialize parameters
  accuracy = dx*relative_accuracy;

  tolerance = 60;
  impatience = 0;

  // initial value, left point and midxle point
  XX2 = x0;
  XX1 = XX2-10*accuracy; // this value 10*accuracy is meanless, just to make sure the check in the while statement goes through

  while (abs(XX2-XX1)>accuracy)
  {
    XX1 = XX2; // copy values

    // value of function at XX
    F0 = (*func)(XX1) - y; // the value of the function at this point

    // decide X1 and X2 for differentiation
    if (XX1>xL+dx)
      X1 = XX1 - dx;
    else
      X1 = xL;

    if (XX1<xR-dx)
      X2 = XX1 + dx;
    else
      X2 = xR;

    // get values at X1 and X2
    F1 = (*func)(X1);
    F2 = (*func)(X2);
    F3 = (F1-F2)/(X1-X2); // derivative at XX1

    XX2 = XX1 - F0/F3; // Newton's mysterious method

    impatience = impatience + 1;
    //cout << "impatience=" << impatience << endl;
    if (impatience>tolerance)
    {
      cout << "invertFunc: " << "max number of iterations reached." << endl;
      exit(-1);
    }

  } // <=> abs(XX2-XX1)>accuracy

  return XX2;
}



//**********************************************************************
vector<double> *zq_x_global, *zq_y_global;
double invertTableDirect_hook(double xx) {return interpCubicDirect(zq_x_global,zq_y_global,xx);}
double invertTableDirect(vector<double>* x, vector<double>* y, double y0, double x0, double relative_accuracy)
// Return x0=y^(-1)(y0) for y=y(x); use interpCubic and invertFunc.
//  -- x,y: the independent and dependent variables. x is assumed to be equal-spaced.
//  -- y0: where the inversion should be performed.
//  -- x0: initial guess
{
  long size = x->size();
  if (size==1) return (*y)[0];
  zq_x_global = x; zq_y_global = y;
  return invertFunc(&invertTableDirect_hook, y0, (*x)[0], (*x)[size-1], (*x)[1]-(*x)[0], x0, relative_accuracy);
}




//**********************************************************************
vector<double> stringToDoubles(string str)
// Return a vector of doubles from the string "str". "str" should
// be a string containing a line of data.
{
  stringstream sst(str+" "); // add a blank at the end so the last data will be read
  vector<double> valueList;
  double val;
  sst >> val;
  while (sst.eof()==false)
  {
    valueList.push_back(val);
    sst >> val;
  }
  return valueList;
}


//**********************************************************************
double stringToDouble(string str)
// Return the 1st doubles number read from the string "str". "str" should be a string containing a line of data.
{
  stringstream sst(str+" "); // add a blank at the end so the last data will be read
  double val;
  sst >> val;
  return val;
}



//**********************************************************************
vector< vector<double>* >* readBlockData(istream &stream_in)
// Return a nested vector of vector<double>* object. Each column of data
// is stored in a vector<double> array and the collection is the returned
// object. Data are read from the input stream "stream_in". Each line
// of data is processed by the stringToDoubles function. Note that the
// data block is dynamicall allocated and is not release within the
// function.
// Note that all "vectors" are "new" so don't forget to delete them.
// Warning that also check if the last line is read correctly. Some files
// are not endded properly and the last line is not read.
{
  vector< vector<double>* >* data;
  vector<double> valuesInEachLine;
  long lineSize;
  long i; // temp variable
  char buffer[99999]; // each line should be shorter than this

  // first line:
  stream_in.getline(buffer,99999);
  valuesInEachLine = stringToDoubles(buffer);
  // see if it is empty:
  lineSize = valuesInEachLine.size();
  if (lineSize==0)
  {
    // empty:
    cout << "readBlockData warning: input stream has empty first row; no data read" << endl;
    return NULL;
  }
  else
  {
    // not empty; allocate memory:
    data = new vector< vector<double>* >(lineSize);
    for (i=0; i<lineSize; i++) (*data)[i] = new vector<double>;
  }

  // rest of the lines:
  while (stream_in.eof()==false)
  {
    // set values:
    for (i=0; i<lineSize; i++) (*(*data)[i]).push_back(valuesInEachLine[i]);
    // next line:
    stream_in.getline(buffer,99999);
    valuesInEachLine = stringToDoubles(buffer);
  }

  return data;
}


//**********************************************************************
void releaseBlockData(vector< vector<double>* >* data)
// Use to delete the data block allocated by readBlockData function.
{
  if (data)
  {
    for (unsigned long i=0; i<data->size(); i++) delete (*data)[i];
    delete data;
  }
}


//**********************************************************************
// From Wikipedia --- the free encyclopeida
//
// Recursive auxiliary function for adaptiveSimpsons() function below
//
double adaptiveSimpsonsAux(double (*f)(double), double a, double b, double epsilon,
                           double S, double fa, double fb, double fc, int bottom) {
  double c = (a + b)/2, h = b - a;
  double d = (a + c)/2, e = (c + b)/2;
  double fd = f(d), fe = f(e);
  double Sleft = (h/12)*(fa + 4*fd + fc);
  double Sright = (h/12)*(fc + 4*fe + fb);
  double S2 = Sleft + Sright;
  if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)
    return S2 + (S2 - S)/15;
  return adaptiveSimpsonsAux(f, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +
         adaptiveSimpsonsAux(f, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);
}
//
// Adaptive Simpson's Rule
//
double adaptiveSimpsons(double (*f)(double),   // ptr to function
                        double a, double b,  // interval [a,b]
                        double epsilon,  // error tolerance
                        int maxRecursionDepth) {   // recursion cap
  double c = (a + b)/2, h = b - a;
  double fa = f(a), fb = f(b), fc = f(c);
  double S = (h/6)*(fa + 4*fc + fb);
  return adaptiveSimpsonsAux(f, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);
}


//**********************************************************************
double qiu_simpsons(double (*f)(double), // ptr to function
                    double a, double b, // interval [a,b]
                    double epsilon, int maxRecursionDepth) // recursion maximum
// My version of the adaptive simpsons integration method.
{
  double f_1=f(a)+f(b), f_2=0., f_4=0.; // sum of values of f(x) that will be weighted by 1, 2, 4 respectively, depending on where x is
  double sum_previous=0., sum_current=0.; // previous and current sum (intgrated value)

  long count = 1; // how many new mid-points are there
  double length = (b-a), // helf length of the interval
         step = length/count; // mid-points are located at a+(i+0.5)*step, i=0..count-1

  int currentRecursionDepth = 1;

  f_4 = f(a+0.5*step); // mid point of [a,b]
  sum_current = (length/6)*(f_1 + f_2*2. + f_4*4.); // get the current sum

  do
  {
    sum_previous = sum_current; // record the old sum
    f_2 += f_4; // old mid-points with weight 4 will be new mid-points with weight 2

    count*=2; // increase number of mid-points
    step/=2.0; // decrease jumping step by half
    f_4 = 0.; // prepare to sum up f_4
    for (int i=0; i<count; i++) f_4 += f(a+step*(i+0.5)); // sum up f_4

    sum_current = (length/6/count)*(f_1 + f_2*2. + f_4*4.); // calculate current sum
    //cout << sum_current << endl;

    if (currentRecursionDepth>maxRecursionDepth)
    {
      cout << endl << "Warning qiu_simpsons: maximum recursion depth reached!" << endl << endl;
      break; // safety treatment
    }
    else currentRecursionDepth++;

  } while (abs(sum_current-sum_previous)>epsilon);

  return sum_current;
}

//**********************************************************************
double qiu_simpsonsRel(double (*f)(double), // ptr to function
                    double a, double b, // interval [a,b]
                    double epsilon, int maxRecursionDepth) // recursion maximum
// My version of the adaptive simpsons integration method.
{
  double f_1=f(a)+f(b), f_2=0., f_4=0.; // sum of values of f(x) that will be weighted by 1, 2, 4 respectively, depending on where x is
  double sum_previous=0., sum_current=0.; // previous and current sum (intgrated value)

  long count = 1; // how many new mid-points are there
  double length = (b-a), // helf length of the interval
         step = length/count; // mid-points are located at a+(i+0.5)*step, i=0..count-1

  int currentRecursionDepth = 1;

  f_4 = f(a+0.5*step); // mid point of [a,b]
  sum_current = (length/6)*(f_1 + f_2*2. + f_4*4.); // get the current sum

  do
  {
    sum_previous = sum_current; // record the old sum
    f_2 += f_4; // old mid-points with weight 4 will be new mid-points with weight 2

    count*=2; // increase number of mid-points
    step/=2.0; // decrease jumping step by half
    f_4 = 0.; // prepare to sum up f_4
    for (int i=0; i<count; i++) f_4 += f(a+step*(i+0.5)); // sum up f_4

    sum_current = (length/6/count)*(f_1 + f_2*2. + f_4*4.); // calculate current sum
    //cout << sum_current << endl;

    if (currentRecursionDepth>maxRecursionDepth)
    {
      cout << endl << "Warning qiu_simpsons: maximum recursion depth reached!" << endl << endl;
      break; // safety treatment
    }
    else currentRecursionDepth++;

  } while (abs(sum_current-sum_previous)/(sum_current-sum_previous)>epsilon);

  return sum_current;
}


//**********************************************************************
string toLower(string str)
// Convert all character in string to lower case
{
  string tmp = str;
  for (string::iterator it=tmp.begin(); it<=tmp.end(); it++) *it = tolower(*it);
  return tmp;
}

//**********************************************************************
string trim(string str)
// Convert all character in string to lower case
{
  string tmp = str;
  long number_of_char = 0;
  for (size_t ii=0; ii<str.size(); ii++)
    if (str[ii]!=' ' && str[ii]!='\t')
    {
      tmp[number_of_char]=str[ii];
      number_of_char++;
    }
  tmp.resize(number_of_char);
  return tmp;
}


//**********************************************************************
long binarySearch(vector<double>* A, double value, bool skip_out_of_range)
// Return the index of the largest number less than value in the list A
// using binary search. Index starts with 0.
// If skip_out_of_range is set to true, then it will return -1 for those
// samples that are out of the table range.
{
   int length = A->size();
   int idx_i, idx_f, idx;
   idx_i = 0;
   idx_f = length-1;
   if(value > (*A)[idx_f])
   {
      if (skip_out_of_range) return -1;
      cout << "binarySearch: desired value is too large, exceeding the end of the table." << endl;
      exit(-1);
   }
   if(value < (*A)[idx_i])
   {
      if (skip_out_of_range) return -1;
      cout << "binarySearch: desired value is too small, exceeding the beginning of table." << endl;
      exit(-1);
   }
   idx = (int) floor((idx_f+idx_i)/2.);
   while((idx_f-idx_i) > 1)
   {
     if((*A)[idx] < value)
        idx_i = idx;
     else
        idx_f = idx;
     idx = (int) floor((idx_f+idx_i)/2.);
   }
   return(idx_i);
}


//**********************************************************************
long double gamma_function(long double x)
// gamma.cpp -- computation of gamma function.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
// Returns gamma function of argument 'x'.
//
// NOTE: Returns 1e308 if argument is a negative integer or 0,
//      or if argument exceeds 171.
{
  int i,k,m;
  long double ga,gr,r=0,z;

  static long double g[] = {
    1.0,
    0.5772156649015329,
   -0.6558780715202538,
   -0.420026350340952e-1,
    0.1665386113822915,
   -0.421977345555443e-1,
   -0.9621971527877e-2,
    0.7218943246663e-2,
   -0.11651675918591e-2,
   -0.2152416741149e-3,
    0.1280502823882e-3,
   -0.201348547807e-4,
   -0.12504934821e-5,
    0.1133027232e-5,
   -0.2056338417e-6,
    0.6116095e-8,
    0.50020075e-8,
   -0.11812746e-8,
    0.1043427e-9,
    0.77823e-11,
   -0.36968e-11,
    0.51e-12,
   -0.206e-13,
   -0.54e-14,
    0.14e-14};

  if (x > 171.0) return 1e308;    // This value is an overflow flag.
  if (x == (int)x) {
    if (x > 0.0) {
      ga = 1.0;               // use factorial
      for (i=2;i<x;i++) {
        ga *= i;
      }
     }
     else
      ga = 1e308;
   }
   else {
    if (fabs(x) > 1.0) {
      z = fabs(x);
      m = (int)z;
      r = 1.0;
      for (k=1;k<=m;k++) {
        r *= (z-k);
      }
      z -= m;
    }
    else
      z = x;
    gr = g[24];
    for (k=23;k>=0;k--) {
      gr = gr*z+g[k];
    }
    ga = 1.0/(gr*z);
    if (fabs(x) > 1.0) {
      ga *= r;
      if (x < 0.0) {
        ga = -M_PI/(x*ga*sin(M_PI*x));
      }
    }
  }
  return ga;
}






//**********************************************************************
long double log_gamma_function(long double x)
// Return ln(Gamma(x)).
// Courtesy to Codecog.org.
// This funcion is under GP license.
{
  if (x <= 0 || x > XBIG)
    return HUGE_VAL;

  if (x <= MACHINE_EPSILON)
    return -log(x);

  if (x <= 4)
  {
    double
    p1[8] =
    {
      4.945235359296727046734888E+00, 2.018112620856775083915565E+02, 2.290838373831346393026739E+03,
      1.131967205903380828685045E+04, 2.855724635671635335736389E+04, 3.848496228443793359990269E+04,
      2.637748787624195437963534E+04, 7.225813979700288197698961E+03
    },
    q1[8] =
    {
      6.748212550303777196073036E+01, 1.113332393857199323513008E+03, 7.738757056935398733233834E+03,
      2.763987074403340708898585E+04, 5.499310206226157329794414E+04, 6.161122180066002127833352E+04,
      3.635127591501940507276287E+04, 8.785536302431013170870835E+03
    },
    p2[8] =
    {
      4.974607845568932035012064E+00, 5.424138599891070494101986E+02, 1.550693864978364947665077E+04,
      1.847932904445632425417223E+05, 1.088204769468828767498470E+06, 3.338152967987029735917223E+06,
      5.106661678927352456275255E+06, 3.074109054850539556250927E+06
    },
    q2[8] =
    {
      1.830328399370592604055942E+02, 7.765049321445005871323047E+03, 1.331903827966074194402448E+05,
      1.136705821321969608938755E+06, 5.267964117437946917577538E+06, 1.346701454311101692290052E+07,
      1.782736530353274213975932E+07, 9.533095591844353613395747E+06
    };

    double corr = x >= PNT68 ? 0 : -log(x);
    double xden = 1, xnum = 0, xm = x <= 1.5 ? (x > 0.5 ? x - 1 : x) : x - 2;
    bool flag = false;
    if (x <= 1.5 && (x <= 0.5 || x >= PNT68)) flag = true;
    double *p = flag ? p1 : p2, *q = flag ? q1 : q2;

    xnum = xnum * xm + p[0], xden = xden * xm + q[0];
    xnum = xnum * xm + p[1], xden = xden * xm + q[1];
    xnum = xnum * xm + p[2], xden = xden * xm + q[2];
    xnum = xnum * xm + p[3], xden = xden * xm + q[3];
    xnum = xnum * xm + p[4], xden = xden * xm + q[4];
    xnum = xnum * xm + p[5], xden = xden * xm + q[5];
    xnum = xnum * xm + p[6], xden = xden * xm + q[6];
    xnum = xnum * xm + p[7], xden = xden * xm + q[7];

    return (x > 1.5 ? 0 : corr) + xm * ((flag ? D1 : D2) + xm * (xnum / xden));
  }

  if (x <= 12)
  {
    double xm = x - 4, xden = -1, xnum = 0,
    p[8] =
    {
      1.474502166059939948905062E+04, 2.426813369486704502836312E+06, 1.214755574045093227939592E+08,
      2.663432449630976949898078E+09, 2.940378956634553899906876E+010,1.702665737765398868392998E+011,
      4.926125793377430887588120E+011, 5.606251856223951465078242E+011
    },
    q[8] =
    {
      2.690530175870899333379843E+03, 6.393885654300092398984238E+05, 4.135599930241388052042842E+07,
      1.120872109616147941376570E+09, 1.488613728678813811542398E+010, 1.016803586272438228077304E+011,
      3.417476345507377132798597E+011, 4.463158187419713286462081E+011
    };

    xnum = xnum * xm + p[0], xden = xden * xm + q[0];
    xnum = xnum * xm + p[1], xden = xden * xm + q[1];
    xnum = xnum * xm + p[2], xden = xden * xm + q[2];
    xnum = xnum * xm + p[3], xden = xden * xm + q[3];
    xnum = xnum * xm + p[4], xden = xden * xm + q[4];
    xnum = xnum * xm + p[5], xden = xden * xm + q[5];
    xnum = xnum * xm + p[6], xden = xden * xm + q[6];
    xnum = xnum * xm + p[7], xden = xden * xm + q[7];

    return D4 + xm * (xnum / xden);
  }

  double res = 0;
  if (x <= FRTBIG)
  {
    double xsq = x * x,
    c[7] =
    {
      -1.910444077728E-03, 8.4171387781295E-04, -5.952379913043012E-04, 7.93650793500350248E-04,
      -2.777777777777681622553E-03, 8.333333333333333331554247E-02, 5.7083835261E-03
    };
    res = c[6];
    res = res / xsq + c[0];
    res = res / xsq + c[1];
    res = res / xsq + c[2];
    res = res / xsq + c[3];
    res = res / xsq + c[4];
    res = res / xsq + c[5];
  }
  double corr = log(x);
  res /= x, res += SQRTPI - corr / 2 + x * (corr - 1);
  return res;
}





//**********************************************************************
double beta_function(double x, double y)
// B(x,y):=Gamma(x)Gamma(y)/Gamma(x+y)
{
    return exp(log_gamma_function(x)+log_gamma_function(y)-log_gamma_function(x+y));
}





//**********************************************************************
double binomial_coefficient(double n, double k)
// Returns the binomial coefficient
// C_n^k := Gamma(n+1) / (Gamma(k+1)*Gamma(n-k+1))
{
    return 1.0/(n+1)/beta_function(k+1,n-k+1);
}





//**********************************************************************
void print_progressbar(double percentage, int length, string symbol)
// Print out a progress bar with the given percentage. Use a negative value to reset the progress bar.
{
  static int status=0;
  static int previous_stop=0;

  if (percentage<0)
  {
    // reset
    status = 0;
    previous_stop = 0;
  }

  // initializing progressbar
  if (status==0)
  {
    cout << "\r";
    cout << "[";
    for (int i=1; i<=length; i++) cout << " ";
    cout << "]";
    cout << "\r";
    cout << "[";
  }

  // plot status
  int stop;
  if (percentage>=0.99) stop=0.99*length;
  else stop = percentage*length;
  for (int i=previous_stop; i<stop; i++) cout << symbol;
  if (previous_stop<stop) previous_stop=stop;

  // simulate a rotating bar
  if (status==0) cout << "-";
  switch (status)
  {
    case 1: cout << "\\"; break;
    case 2: cout << "|"; break;
    case 3: cout << "/"; break;
    case 4: cout << "-"; break;
  }
  cout << "\b";
  status++;
  if (status==5) status=1;
  cout.flush();
}


//**********************************************************************
void formatedPrint(ostream& os, int count, ...)
// For easier scientific data outputing.
{
  va_list ap;
  va_start(ap, count); //Requires the last fixed parameter (to get the address)
  for(int j=0; j<count; j++)
      os << scientific << setprecision(OUTPUT_PRECISION) << "  " << va_arg(ap, double); //Requires the type to cast to. Increments ap to the next argument.
  va_end(ap);
  os << endl;
}



//**********************************************************************
void display_logo(int which)
// Personal amusement.
{
  switch (which)
  {
    case 1:
    cout << " ____  ____            _                    " << endl;
    cout << "|_   ||   _|          (_)                   " << endl;
    cout << "  | |__| |    .---.   __    _ .--.    ____  " << endl;
    cout << "  |  __  |   / /__\\\\ [  |  [ `.-. |  [_   ] " << endl;
    cout << " _| |  | |_  | \\__.,  | |   | | | |   .' /_ " << endl;
    cout << "|____||____|  '.__.' [___] [___||__] [_____]" << endl;
    cout << "                                            " << endl;
    break;

    case 2:
    cout << ":::    ::: :::::::::: ::::::::::: ::::    ::: :::::::::" << endl;
    cout << ":+:    :+: :+:            :+:     :+:+:   :+:      :+: " << endl;
    cout << "+:+    +:+ +:+            +:+     :+:+:+  +:+     +:+  " << endl;
    cout << "+#++:++#++ +#++:++#       +#+     +#+ +:+ +#+    +#+   " << endl;
    cout << "+#+    +#+ +#+            +#+     +#+  +#+#+#   +#+    " << endl;
    cout << "#+#    #+# #+#            #+#     #+#   #+#+#  #+#     " << endl;
    cout << "###    ### ########## ########### ###    #### #########" << endl;
    break;

    case 3:
    cout << " __  __     ______     __     __   __     _____    " << endl;
    cout << "/\\ \\_\\ \\   /\\  ___\\   /\\ \\   /\\ '-.\\ \\   /\\___  \\  " << endl;
    cout << "\\ \\  __ \\  \\ \\  __\\   \\ \\ \\  \\ \\ \\-.  \\  \\/_/  /__ " << endl;
    cout << " \\ \\_\\ \\_\\  \\ \\_____\\  \\ \\_\\  \\ \\_\\\\'\\_\\   /\\_____\\" << endl;
    cout << "  \\/_/\\/_/   \\/_____/   \\/_/   \\/_/ \\/_/   \\/_____/" << endl;
    break;

  }

}



//**********************************************************************
void GaussLegendre_getWeight(int npts,double* xg,double* wg, double A, double B, int iop)
// Calculate the sampling location and weight for Gauss-Legendre quadrature
// -- From Hirano and Nara's MC-KLN code.
{
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//  gauss.f: Points and weights for Gaussian quadrature                 c
//                       c
//  taken from: "Projects in Computational Physics" by Landau and Paez  c
//         copyrighted by John Wiley and Sons, New York            c
//                                                                      c
//  written by: Oregon State University Nuclear Theory Group            c
//         Guangliang He & Rubin H. Landau                         c
//  supported by: US National Science Foundation, Northwest Alliance    c
//                for Computational Science and Engineering (NACSE),    c
//                US Department of Energy                          c
//                       c
//  comment: error message occurs if subroutine called without a main   c
//  comment: this file has to reside in the same directory as integ.c   c
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  static const double EPS = 3.0e-14;
  int m=(npts+1)/2;
  for(int i=0;i<m;i++) {
    double  t=cos(M_PI*(i+1-0.25)/(npts+0.5));
    double t1=t;
    double pp;
    do {
        double p1=1.0;
        double p2=0.0;
        double aj=0.0;
        for(int j=0;j<npts;j++) {
            double p3=p2;
            p2=p1;
            aj += 1.0;
            p1=((2.0*aj-1.0)*t*p2-(aj-1.0)*p3)/aj;
        }
        pp=npts*(t*p1-p2)/(t*t-1.0);
        t1=t;
        t=t1-p1/pp;
    } while(abs(t-t1)>EPS);
    xg[i]=-t;
    xg[npts-1-i]=t;
    wg[i]=2.0/((1.0-t*t)*pp*pp);
    wg[npts-1-i]=wg[i];
    }

//GaussLegendre::GaussRange(int N,int iop,double A,double B,
//  double* xg1,double* wg1)
//{
//     transform gausspoints to other range than [-1;1]
//     iop = 1  [A,B]       uniform
//     iop = 2  [0,inf]     A is midpoint
//     opt = 3  [-inf,inf]  scale is A
//     opt = 4  [B,inf]     A+2B is midoint
//     opt = 5  [0,B]     AB/(A+B)+ is midoint

  int N=npts;
  double xp, wp;
  for(int i=0; i<N; i++) {
      if(iop == 1) {
          //...... A to B
          xp=(B+A)/2+(B-A)/2*xg[i];
          wp=(B-A)/2*wg[i];
      } else if(iop == -1) {
          //......   A to B
          xp=(B+A)/2+(B-A)/2*xg[i];
          if(i <= N/2)
            xp=(A+B)/2-(xp-A);
          else
            xp=(A+B)/2+(B-xp);
          wp=(B-A)/2*wg[i];
      } else if(iop == 2) {
          //...... zero to infinity
          xp=A*(1+xg[i])/(1-xg[i]);
          double tmp=(1-xg[i]);
          wp=2.*A/(tmp*tmp)*wg[i];
      } else if(iop ==  3) {
          //...... -inf to inf scale A
          xp=A*(xg[i])/(1-xg[i]*xg[i]);
          double tmp=1-xg[i]*xg[i];
          wp=A*(1+xg[i]*xg[i])/(tmp*tmp)*wg[i];
      } else if(iop == 4) {
          //......  B to inf,  A+2B is midoint
          xp=(A+2*B+A*xg[i])/(1-xg[i]);
          double tmp=1-xg[i];
          wp=2.*(B+A)/(tmp*tmp)*wg[i];
      } else if(iop == 5) {
          //...... -A to A , scale B
          //xp=A*pow(abs(xg[i]),B) *sign(1.0,xg(i));
          double tmp = xg[i] >= 0 ? 1.0 : -1.0;
          xp=A*pow(abs(xg[i]),B) * tmp;
          //xp=A*pow(abs(xg[i]),B) *sign(1.0,xg(i));
          wp=A*B*pow(abs(xg[i]),(B-1))*wg[i];
      } else if(iop ==  6) {
          //...... 0 to B , AB/(A+B) is midpoint
          xp=A*B*(1+xg[i])/(B+A-(B-A)*xg[i]);
          double tmp = B+A-(B-A)*xg[i];
          wp=2*A*B*B/(tmp*tmp)*wg[i];
      } else {
          cerr << " invalid option iop = " << iop << endl;
          exit(-1);
      }
      xg[i]=xp;
      wg[i]=wp;
  }
}





//**********************************************************************
void get_bin_average_and_count(istream& is, ostream& os, vector<double>* bins, long col_to_bin, void (*func)(vector<double>*), long wanted_data_columns, bool silence)
// Group data into bins by set by the "bins". The data in the column
// "col_to_bin" read from "is" are the ones used to determine the binning.
// Once the binning is decided, the averages of all data are calculated.
// The result, which has the same number of rows and the number of bins,
// will be outputted to "os". The i-th line of the output data contains the
// average of data from each column, for the i-th bin. The output will
// have 2 more columns; the 1st being the number count N and the 2nd being
// dN/dx where dx is the bin width.
// The values given in "bins" define the boundaries of bins and is assumed
// to be non-decreasing.
// After each line is decided to go into which bin, the function specified
// by "func" will be called to transform data. It is the transformed data
// that will be averaged. The transformed data can have different number of
// columns than the data passed in, in which case its number of columns
// is specified by "wanted_data_columns". The counting info will still
// be always written in the last two columns.
// The function specified by "func" should accepts a vector of doubles
// which is one line of data, and then modify it as returned result. The
// data vector passed in already has the desired length so it can be modified
// directly.
// The argument "col_to_bin" starts with 1.
// Refer also to getBinnedAverageAndCount MATLAB program.
{
  // initialization
  char* buffer = new char[99999]; // each line should be shorter than this
  // get first line and continue initialization
  is.getline(buffer, 99999);
  vector<double> line_data = stringToDoubles(buffer);
  long number_of_cols = line_data.size();
  long number_of_bins = bins->size()-1;

  if (number_of_cols==0)
  {
    cout << "get_bin_average_and_count error: the input data is empty!" << endl;
    exit(-1);
  }

  // create the counting array
  if (wanted_data_columns>0) number_of_cols = wanted_data_columns;
  double bin_total_and_count[number_of_bins][number_of_cols+2];
  for (long i=0; i<number_of_bins; i++)
  for (long j=0; j<number_of_cols+2; j++)
  {
    bin_total_and_count[i][j] = 0;
  }

  // add up all data
  long number_of_lines=1;
  while (is.eof()==false)
  {
    // determine which bin
    long bin_idx = binarySearch(bins, line_data[col_to_bin-1], true);
    if (bin_idx!=-1)
    {
      // transform data
      line_data.resize(number_of_cols, 0);
      if (func) (*func) (&line_data);
      // add to the counting matrix
      for (long j=0; j<number_of_cols; j++) bin_total_and_count[bin_idx][j] += line_data[j];
      // also the counting column
      bin_total_and_count[bin_idx][number_of_cols] ++;
    }
    // next iteration
    is.getline(buffer, 99999);
    line_data = stringToDoubles(buffer);
    if (number_of_lines % 100000 == 0 && !silence) cout << "Line " << number_of_lines << " reached." << endl;
    number_of_lines++;
  }

  // find the averages
  for (long i=0; i<number_of_bins; i++)
  for (long j=0; j<number_of_cols; j++)
  {
    if (bin_total_and_count[i][number_of_cols]<1e-15) continue;
    bin_total_and_count[i][j] /= bin_total_and_count[i][number_of_cols];
  }


  // get dN/(d bin_width)
  for (long i=0; i<number_of_bins; i++)
  {
    if (bin_total_and_count[i][number_of_cols]<1e-15) continue;
    bin_total_and_count[i][number_of_cols+1] = bin_total_and_count[i][number_of_cols]/((*bins)[i+1]-(*bins)[i]);
  }

  // output
  for (long i=0; i<number_of_bins; i++)
  {
    for (long j=0; j<number_of_cols+2; j++)
    {
      os << scientific << scientific << setprecision(OUTPUT_PRECISION) << bin_total_and_count[i][j] << "  ";
    }
    os << endl;
  }

}
