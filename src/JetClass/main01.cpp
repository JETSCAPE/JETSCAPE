//
//  main01.cpp
//  
//
//  Created by Abhijit Majumder on 11/25/16.
//  A simple program to create a jet using a predetermined
//  algorithm for parton propagation and splitting. 

#include <iostream>
#include <complex>
#include <fstream>
#include <math.h>
#include <assert.h>
#include "JetClass.hpp"


using namespace std;

int main()
{
    double virt;
    double pAssign[4], xLoc[4];
    int i;
    
    double deltaTimeStep;
    
    deltaTimeStep = 0.1;
    
    for (i=0;i<=3; i++) {
        xLoc[i] = 0.0;
    };
    
    pAssign[0] = 11.0;
    
    pAssign[3] = 10.0;
    
    pAssign[1] = pAssign[2] = 1.0;
    
    FourVector p_in(pAssign);
    
    
    cout << " p_in = " << p_in.t() << "  " << p_in.x() << "  " << p_in.y() << "  " << p_in.z() << endl;
 
     Jet jet(p_in);
    
 //   cout << " virt = " << virt << endl ;
    
    cout << " jet eta = " << jet.get_jet_eta() << endl;
    
    
    
    
    
}

/* TDOD
virtual double parton::generate_t()
{
    
    double s_approx = 0.01
    double pi = 3.1415
    
    t = 2.0;
    
    return (t);
    
    
}

*/



