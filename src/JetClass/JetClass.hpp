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
#include "four_vector.hpp"


class Jet
{
    friend class Parton;
    
public:
    
    Jet (double p_in[4]);
    
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


class Parton
{

public:
    
    Parton (int label, int id, int stat, double p[4], double x[4]);
    
    void set_label(int label)
    {
        plabel_ = label;
    };

    void set_id(int id)
    {
        pid_ = id;
    };
    
    void set_stat(int stat)
    {
        pstat_ = stat;
    };
    
    void set_mass(double mass_input)
    {
        mass_ = mass_input;
    }
    
    void set_p(double p[4])
    {
        FourVector p_in_(p);
    };
    
    void set_x(double x[4])
    {
        FourVector x_in_(x);
    };
    
    
    double mass()
    {
        return(mass_);
    }

    
    double get_p(int i)
    {
        return (p_in_.comp(i));
    };
    
    double pl()
    {
        return(std::sqrt( p_in_.x()*p_in_.x() + p_in_.y()*p_in_.y() + p_in_.z()*p_in_.z() ) );
    };

    
    void set_mean_form_time ()
    {
        mean_form_time_ = this->pl()/2/t_;
    }

    
    double get_t ()
    {
        
        t_ = generate_t();
        
        return (t_) ;
    }; // virtuality of particle
    
    virtual double generate_t()
    {
        return (1);
    };
    
    
private:
    
    int pid_;  // particle id ()
    int pstat_; // status of particle
    int plabel_; // the line number in the event record
    FourVector p_in_, x_in_; // internal momentum and position of particle
    double Energy_, t_; // Energy, and t is the standard virtuality variable
    double mean_form_time_ ;  // Mean formation time
    double form_time_ ; //event by event formation time
    double mass_; //mass of the parton

};





#endif /* JetClass_hpp */
