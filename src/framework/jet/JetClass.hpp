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
#include "fjcore.hh"

#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>

class Parton;
class VertexBase;
class FourVector;


/**************************************************************************************************/

//  JET CLASS

/*************************************************************************************************/

//dummy for now figure out after graph strcuture ...

class Jet
{
  Jet() {};
  ~Jet() {};
};


/**************************************************************************************************/

//  PARTON CLASS

/*************************************************************************************************/

class Parton : public fjcore::PseudoJet
{

public:

  Parton() : PseudoJet() {};
  Parton (int label, int id, int stat, double p[4], double x[4]);
  Parton (int label, int id, int stat, double pt, double eta, double phi, double e);
  Parton (int label, int id, int stat, double pt, double eta, double phi, double e, double x[4]);
	  
  virtual ~Parton();
  
    void clear()
    {
        plabel_ = 0;
        pid_ = 0;
        pstat_ = 0;        
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

  // from PseudoJet ... (well keep for now);
    void set_mass(double mass_input)
    {
        mass_ = mass_input;
    };
    
    void set_p(double p[4])
    {
      //FourVector p_in_(p); // error: creates new vector and hence not accessible via class p_in_
      reset_momentum(p[0],p[1],p[2],p[3]);
    };

  // not needed in graph structure
  
    void set_x(double x[4])
    {
      //FourVector
	  x_in_.Set(x);
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

  /*
    int pparent_label()
    {
        return(pparent_label_);
    }
    */
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

  // just operator of PseudoJet ...
  /*
    double get_p(int i)
    {
        return (p_in_.comp(i));
    };
  */
  
    double pl()
    {
      return 0;//(std::sqrt( p_in_.x()*p_in_.x() + p_in_.y()*p_in_.y() + p_in_.z()*p_in_.z() ) );
    }
    
    double get_t ()
    {
        
        t_ = generate_t();
        
        return (t_);
    }; // virtuality of particle
    
    virtual double generate_t()
    {
        return (1);
    };
    
    Parton& operator=(Parton &c)
    {
      //FourVector x_in_;
      
      pid_ = c.pid() ;
      pstat_ = c.pstat() ;
      plabel_ = c.plabel() ;
      //pparent_label_ = c.pparent_label();
      //p_in_ = c.p_in() ;
      //Parton::set_x(c.x_in());

      x_in_ = c.x_in() ;
      
      mass_ = c.mass();
      t_ = c.get_t();
      
      return *this;
    }

   Parton& operator=(const Parton &c)
    {
      //FourVector x_in_;
      
      pid_ = c.pid_;
      pstat_ = c.pstat_ ;
      plabel_ = c.plabel_;
      //pparent_label_ = c.pparent_label();
      //p_in_ = c.p_in() ;
      //Parton::set_x(c.x_in());

      x_in_ = c.x_in_;
      
      mass_ = c.mass_;
      t_ = c.t_;
      
      return *this;
    }
  /*
  friend ostream Print()
  {
    ostream output;
    return output;
  }
  */
  
  /*
  friend ostream &operator<<( ostream &output, 
         Parton & parton ) {
    output<<" Parton: ";
    output << " ID : "<<parton.pid()<<endl;
    //output << "vec(p) = "<<parton.p_in().x()<<" "<<parton.p_in().y()<<" "<<parton.p_in().z()<<" "<<parton.x_in().t();
    //output << " vec(p) = "<<parton.get_p(0)<<" "<<parton.get_p(1)<<" "<<parton.get_p(2)<<" "<<parton.get_p(3)<<endl;
    //output << " vec(px,py,px,e)      = "<<parton(0)<<" "<<parton(1)<<" "<<parton(2)<<" "<<parton(3)<<endl;
    output << " vec(pT,eta,phi,e)    = "<<parton.pt()<<" "<<parton.rap()<<" "<<parton.phi()<<" "<<parton.e()<<endl;
    output << " vec(x,y,z,t)         = "<<parton.x_in().x()<<" "<<parton.x_in().y()<<" "<<parton.x_in().z()<<" "<<parton.x_in().t()<<endl;
    return output;            
      }
  */

  friend ostream &operator<<( ostream &output, 
         Parton & parton ) {
    output<<parton.plabel()<<" "<<parton.pid()<<" "<<parton.pstat()<<" ";
    output<<parton.pt()<<" "<<parton.rap()<<" "<<parton.phi()<<" "<<parton.e()<<" ";
    output<<parton.x_in().x()<<" "<<parton.x_in().y()<<" "<<parton.x_in().z()<<" "<<parton.x_in().t();//<<endl;
    
    //output << "vec(p) = "<<parton.p_in().x()<<" "<<parton.p_in().y()<<" "<<parton.p_in().z()<<" "<<parton.x_in().t();
    //output << " vec(p) = "<<parton.get_p(0)<<" "<<parton.get_p(1)<<" "<<parton.get_p(2)<<" "<<parton.get_p(3)<<endl;
    //output << " vec(px,py,px,e)      = "<<parton(0)<<" "<<parton(1)<<" "<<parton(2)<<" "<<parton(3)<<endl;
    //output << " vec(pT,eta,phi,e)    = "<<parton.pt()<<" "<<parton.rap()<<" "<<parton.phi()<<" "<<parton.e()<<endl;
    //output << " vec(x,y,z,t)         = "<<parton.x_in().x()<<" "<<parton.x_in().y()<<" "<<parton.x_in().z()<<" "<<parton.x_in().t()<<endl;
    return output;            
      }
  
private:
  
  int pid_                ; // particle id ()
  int pstat_              ; // status of particle
  int plabel_             ; // the line number in the event record
  double t_;
  double mean_form_time_  ; // Mean formation time
  double form_time_       ; //event by event formation time
  double mass_            ; //mass of the parton

  FourVector p_in_;
  FourVector x_in_; // position of particle
  
  // following will be in graph strucure
  //int pparent_label_      ; // line number of parent
  //FourVector p_in_, x_in_ ; // internal momentum and position of particle
  //double Energy_, t_      ; // Energy, and t is the standard virtuality variable  
  //Vertex *start_, *end_    ; // the vertex from where it came and the vertex where it will end.
};






/**************************************************************************************************/

//  VERTEX CLASS

/*************************************************************************************************/

class VertexBase
{
    
public:

  VertexBase() {x_in_.Set(0,0,0,0);}
  VertexBase(double x, double y, double z, double t) {x_in_.Set(x,y,z,t);}
  VertexBase(FourVector &x) {set_location(x);}  
  virtual ~VertexBase();
  
  void set_location(FourVector &x)
    {
      x_in_ = x;
    }
  
  FourVector &x_in()
  {
    return(x_in_);
  }

  friend ostream &operator<<( ostream &output, 
         VertexBase & vertex ) {

    output<<vertex.x_in().x()<<" "<<vertex.x_in().y()<<" "<<vertex.x_in().z()<<" "<<vertex.x_in().t();//<<endl;
    
    return output;
  }
private:
  
  FourVector x_in_        ; //location of the vertex
  // parents and siblings from Graph structure later ...
  
};

// include graph edges/node propbably no need to inherit from VertexBase ...
class Vertex : public VertexBase
{
  
public:
  
  Vertex() : VertexBase() {};
  Vertex(FourVector &x) : VertexBase(x) {};
  virtual ~Vertex();
  
};

#endif /* JetClass_hpp */
