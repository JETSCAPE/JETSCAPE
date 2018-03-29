//
//  JetScapeConstants.h
//  
//
//  Created by Abhijit Majumder on 11/26/16.
//
//

#ifndef JETSCAPECONSTANTS_H
#define JETSCAPECONSTANTS_H

namespace Jetscape {


// define the usual constants
static double pi = 3.141592653589793 ;

static double nf = 3.0;

static double Cf = 4.0/3.0 ;

static double Tf = 0.5 ;

static double Ca = 3.0 ;

static double Nc = 3.0;

static double Lambda_QCD = 0.2;
// 0.4 is the value chosen in JETSET
    
static double fmToGeVinv = 5.0;
/// < should be 1/0.197, but 5 helps in debugging. 

static double zeta3 = 1.20206;

static double mu = 0.722;

/* When the code becomes really accurate, a more accurate value for this can be used  */

/*  the following is the maximum value from the standard C++ random number generator */
static double maxN = double(pow(2.0,31.0) - 1.0) ;

static double a_very_large_number = maxN ;

/* the following 2 lines control the error in the analytical part of the calculation. */
/* Note analytical approximation, cannot be rectified by more statistics             */
/* However, more accurate analytical calculation will require less statistics to obtain smooth distributions */
/* the value is something for the user to choose based on his/her computing resources  */

static double error = 0.02 ;

static double approx = 0.02 ;

static double s_error = 0.001;
static double s_approx = 0.001;

static double E_minimum = 1.0 ;

static double rounding_error = 1e-6; // slightly more than float precision

/**************************************************************************************/

/* the standard PDG particle id codes for the gluon and the d quark */
static int gid=21;

static int qid=1;

static int uid=2;

static int did=1;

static int sid=3;
/*******************************************************************/

};
#endif // JETSCAPECONSTANTS_H
