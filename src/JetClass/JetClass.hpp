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
#include <vector>


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
            jet_p_[i] = p_in[i];
        };
    };
    
    double get_jet_p ()
    {
        return ( sqrt( jet_p_[1]*jet_p_[1] + jet_p_[2]*jet_p_[2] + jet_p_[3]*jet_p_[3] ) );
    }
    
    double get_jet_eta ()
    {
        double p_mod;
        
        p_mod = get_jet_p ();
        
        jet_eta_ = log( ( p_mod + jet_p_[3] )/( p_mod - jet_p_[3] )  )/2.0;
        
        return (jet_eta_);
    }
    
    double get_jet_phi()
    {
        
        if (jet_p_[1]!=0.0)
        {
            jet_phi_ = atan(jet_p_[2]/jet_p_[1]);
        }
        else
        {
            if (jet_p_[2]>0.0)
            {
                jet_phi_ = pi/2.0;
            }
            else
            {
                jet_phi_ = 3.0*pi/2.0;
            }
        }
        return(jet_phi_);
    }
    
private:
    
    double jet_p_[4]; // momenta of jet
    double jet_eta_; //eta of jet
    double jet_phi_; //phi of jet
    
};


class Parton
{

public:
    Parton (int label, int id, int stat, int parent_label, double p[4], double x[4]);
    
    Parton ()
    :plabel_(0)
    ,pid_(0)
    ,pstat_(0)
    {
        p_in_.clear();
        x_in_.clear();
    }
    
    void clear()
    {
        plabel_ = 0;
        pid_ = 0;
        pstat_ = 0;
        pparent_label_ = 0;
        p_in_.clear();
        x_in_.clear();
    }
    
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
    
    void set_parent_label(int parent_label)
    {
        pparent_label_ = parent_label;
    };
    
    void set_mass(double mass_input)
    {
        mass_ = mass_input;
    };
    
    void set_p(double p[4])
    {
        FourVector p_in_(p);
    };
    
    void set_x(double x[4])
    {
        FourVector x_in_(x);
    };
    
    void set_mean_form_time ()
    {
        mean_form_time_ = this->pl()/2/t_;
    };

    
    int pid()
    {
        return(pid_);
    };
    
    int pstat()
    {
        return(pstat_);
    };
    
    int plabel()
    {
        return(plabel_);
    };
    
    int pparent_label()
    {
        return(pparent_label_);
    }
    
    FourVector &p_in()
    {
        return(p_in_);
    }
    
    FourVector &x_in()
    {
        return(x_in_);
    }
    
    double mass()
    {
        return(mass_);
    }

    double mean_form_time()
    {
        return(mean_form_time_);
    }
    
    double get_p(int i)
    {
        return (p_in_.comp(i));
    };
    
    double pl()
    {
        return(std::sqrt( p_in_.x()*p_in_.x() + p_in_.y()*p_in_.y() + p_in_.z()*p_in_.z() ) );
    };
    
    double get_t ()
    {
        
        t_ = generate_t();
        
        return (t_);
    }; // virtuality of particle
    
    virtual double generate_t()
    {
        return (1);
    };
    
    Parton &operator=(Parton &c)
    {
        pid_ = c.pid() ;
        pstat_ = c.pstat() ;
        plabel_ = c.plabel() ;
        pparent_label_ = c.pparent_label();
        p_in_ = c.p_in() ;
        x_in_ = c.x_in() ;
        mass_ = c.mass();
        t_ = c.get_t();
        
        return(*this);
    }
    
private:
    
    int pid_                ; // particle id ()
    int pstat_              ; // status of particle
    int plabel_             ; // the line number in the event record
    int pparent_label_      ; // line number of parent
    FourVector p_in_, x_in_ ; // internal momentum and position of particle
    double Energy_, t_      ; // Energy, and t is the standard virtuality variable
    double mean_form_time_  ; // Mean formation time
    double form_time_       ; //event by event formation time
    double mass_            ; //mass of the parton

};



class Vertex
{
    
    
public:
    Vertex(int num_parents, int num_siblings, FourVector &x, vector<Parton> &parent, vector<Parton> &sibling);
    
    Vertex()
    :Nparents_(0)
    ,Nsiblings_(0)
    {
        x_in_.clear();
        for (int i=0; i<Nparents_; i++)
        {
            parent_[i].clear();
        }
        for (int j=0; j<Nsiblings_; j++)
        {
            sibling_[j].clear();
        }
    };
    
    void set_num_parents(int num_parents)
    {
        Nparents_ = num_parents;
    }
    
    void set_num_siblings(int num_siblings)
    {
        Nsiblings_ = num_siblings ;
    }
    
    void set_location(FourVector &x)
    {
        x_in_ = x;
    }
    
    void set_parents(vector<Parton> &parent, vector<Parton> &sibling)
    {
        for (int i=0; i<Nparents_; i++)
        {
            parent_.push_back(parent[i]);
        }
        for (int j=0; j<Nsiblings_;j++)
        {
            sibling_.push_back(sibling[j]);
        }
    }
    
private:
    int Nparents_, Nsiblings_ ; // number of parents and siblings.
    FourVector x_in_        ; //location of the vertex
    vector<Parton>  parent_, sibling_ ; // partons that connect to a vertex
    
};


#endif /* JetClass_hpp */
