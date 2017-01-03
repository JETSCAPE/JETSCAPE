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


class Jet
{
    friend class parton;
    
public:
    
    Jet (double p_in[4]);
    
    void Set_Jet_p(double p_in[4])
    {
        int i;
        for (i=0; i<=3; i++)
        {
            jet_p[i] = p_in[i];
        };
    };
    
    double Get_Jet_p ()
    {
        return ( sqrt( jet_p[1]*jet_p[1] + jet_p[2]*jet_p[2] + jet_p[3]*jet_p[3] ) );
    }
    
    double Get_Jet_eta ()
    {
        double p_mod;
        
        p_mod = Get_Jet_p ();
        
        jet_eta = log( ( p_mod + jet_p[3] )/( p_mod - jet_p[3] )  )/2.0;
        
        return (jet_eta);
    }
    
    double Get_Jet_phi()
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


class Parton
{

public:
    
    Parton (int label, int id, int stat, double p_in[4], double x[4]);
    
    void Set_label(int label)
    {
        plabel = label;
    };

    void Set_id(int id)
    {
        pid = id;
    };
    
    void Set_stat(int stat)
    {
        pstat = stat;
    };
    
    void Set_p(double p_in[4])
    {
        p[1] = p_in[1];
        p[2] = p_in[2];
        p[3] = p_in[3];
        p[0] = p_in[0];
    };
    
    void Set_x(double x[4])
    {
        int i;
        for (i=0; i<=3 ; i++) {
            xp[i] = x[i];
        };
    };
    
    double Get_p(int i)
    {
        return (p[i]);

    };
    

    void Set_pl()
    {
        pl = std::sqrt( p[1]*p[1] + p[2]*p[2] + p[3]*p[3] );
    };

    
    void Set_m_form_time ()
    {
        m_form_time = pl/2/t;
    }
    
    double Get_t ()
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

class Vertex
{
    
private:
    
    double x[4];
    Parton parent;
    Parton offspring1;
    Parton offspring2;
};



#endif /* JetClass_hpp */
