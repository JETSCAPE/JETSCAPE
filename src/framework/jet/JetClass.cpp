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
#include "constants.h"

Jet::Jet(FourVector p)
{
    set_jet_p(p);
}


Parton::Parton (int label, int id, int stat, int parent_label, double p[4], double x[4])
{
    
    set_label(label);
    
    set_id(id);
    
    set_mass(-1.0);
    switch (id) {
        case 1:  //down quark
        case -1: // anti-down quark
            set_mass(0.01);
            break;
            
        case 2:  // up quark
        case -2:  // anti-up quark
            set_mass(0.005);
            break;
        
        case 3:   // strange quark
        case -3:  // anti-strange quark
            set_mass(0.15);
            break;
            
        case 21: // gluon
            set_mass(0.0);
            break;
            
        default:
        {
            std::cout << " error in id = " << id << std::endl;
            assert(mass()>=0.0);
            break;
        }
    }
    
    set_stat(stat);
    
    set_parent_label(parent_label);

    //DEBUG:
    //cout<<p[0]<<" "<<p[1]<<endl;
    
    set_p(p);

    //DEBUG:
    //cout<<&p_in()<<endl;
    //cout<<get_p(0)<<" "<<get_p(1)<<endl;
    //cout<<p_in().x()<<endl;
    
    set_x(x);
    
    Energy_ = std::sqrt(this->pl()*this->pl() + this->mass()*this->mass());
    

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

