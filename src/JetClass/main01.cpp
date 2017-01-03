//
//  main01.cpp
//  
//
//  Created by Abhijit Majumder on 11/25/16.
//
//

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
    double passign[4], xloc[4];
    int i;
    
    
    for (i=0;i<=3; i++) {
        passign[i] = 1.0;
        xloc[i] = 0.0;
    };
    
    //std::vector<parton> shower;
    Parton quark(1,1,1,passign,xloc);
    
    Jet parentjet(passign);
    
    virt = quark.generate_t();
    
    cout << " virt = " << virt << endl ;
    
    
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



