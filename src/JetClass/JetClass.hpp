//
//  JetClass.hpp
//  
//
//  Created by Abhijit Majumder on 10/6/16.
//
//

#ifndef JetClass_hpp
#define JetClass_hpp

#include <stdio.h>
#include <math.h>
#include "constants.h"


class jet
{
    friend class parton;
    
public:
    
    jet (double p_in[4]);
    
    void set_jet_p(double p_in[4])
    {
        int i;
        for (i=0; i<=3; i++)
        {
            jet_p[i] = p_in[i];
        };
    };
    
    double get_jet_p ()
    {
        return ( sqrt( jet_p[1]*jet_p[1] + jet_p[2]*jet_p[2] + jet_p[3]*jet_p[3] ) );
    }
    
    double get_jet_eta ()
    {
        double p_mod;
        
        p_mod = get_jet_p ();
        
        jet_eta = log( ( p_mod + jet_p[3] )/( p_mod - jet_p[3] )  )/2.0;
        
        return (jet_eta);
    }
    
    double get_jet_phi()
    {
        
        if (jet_p[1]!=0.0)
        {
            jet_phi = atan(jet_p[2]/jet_p[1]);
        }
        else
        {
            if (jet_p[2]>0.0)
            {
                jet_phi = pi/2.0;
            }
            else
            {
                jet_phi = 3.0*pi/2.0;
            }
        }
        return(jet_phi);
    }
    
private:
    
    double jet_p[4]; // momenta of jet
    double jet_eta; //eta of jet
    double jet_phi; //phi of jet
    
};


class parton
{

public:
    
    parton (int label, int id, int stat, double p_in[4], double x[4]);
    
    void set_label(int label)
    {
        plabel = label;
    };

    void set_id(int id)
    {
        pid = id;
    };
    
    void set_stat(int stat)
    {
        pstat = stat;
    };
    
    void set_p(double p_in[4])
    {
        p[1] = p_in[1];
        p[2] = p_in[2];
        p[3] = p_in[3];
        p[0] = p_in[0];
    };
    
    void set_x(double x[4])
    {
        int i;
        for (i=0; i<=3 ; i++) {
            xp[i] = x[i];
        };
    };
    
    double get_p(int i)
    {
        return (p[i]);

    };
    

    void set_pl()
    {
        pl = std::sqrt( p[1]*p[1] + p[2]*p[2] + p[3]*p[3] );
    };

    
    void set_m_form_time ()
    {
        m_form_time = pl/2/t;
    }
    
    double get_t ()
    {
        
        t = generate_t();
        
        return (t) ;
    }; // virtuality of particle
    
    virtual double generate_t()
    {
        return (1);
    };
    
    
private:
    
    int pid;  // particle id ()
    int pstat; // status of particle
    int plabel; // the line number in the event record
    double mass;
    double eng; // energy of particle
    double pj[4]; // momenta of particle in jet frame
    double pp, pm, t ; //light cone momenta
    double p[4]; // momenta of particle in CM frame
    double xp[4];//location of particle in CM frame
    double pl; // modulus of 3-momentum
    double m_form_time ;  // Mean formation time
    double form_time ; //event by event formation time

};





#endif /* JetClass_hpp */
