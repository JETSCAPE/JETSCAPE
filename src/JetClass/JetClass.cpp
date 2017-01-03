//
//  JetClass.cpp
//  
//
//  Created by Abhijit Majumder on 10/6/16.
//
//

#include <iostream>
#include <complex>
#include <fstream>
#include <math.h>
#include <assert.h>
#include "JetClass.hpp"
//#include "constants.h"

Jet::Jet(double p_in[4])
{
    Set_Jet_p(p_in);
}

Parton::Parton (int label, int id, int stat, double p_in[4], double x[4])
{
    
    Set_label(label);
    
    Set_id(id);
    
    mass = -1.0;
    switch (id) {
        case 1:  //down quark
        case -1: // anti-down quark
            mass = 0.01;
            break;
            
        case 2:  // up quark
        case -2:  // anti-up quark
            mass = 0.005;
            break;
        
        case 3:   // strange quark
        case -3:  // anti-strange quark
            mass = 0.15;
            break;
            
        case 0: // gluon
            mass = 0.0;
            break;
            
        default:
        {
            std::cout << " error in id = " << id << std::endl;
            assert(mass>=0.0);
            break;
        }
    }
    
    Set_stat(stat);
    
    Set_p(p_in);
    
    Set_x(x);
    
    Set_pl();
    
    p[0] = std::sqrt(pl*pl + mass*mass);
    
    pp = (p[0] + pl)/std::sqrt(2);
    
    pm = (p[0] - pl)/std::sqrt(2);
    
    
   /* @TODO 
    
    void parton::set_pperpx(jet parentjet)
    {
        
        pperpx = 0.0;
        
        
    };
    
    
    void parton::set_pperpy(jet parentjet)
    {
        pperpy = 0.0;
    };

    */
    
    
}

