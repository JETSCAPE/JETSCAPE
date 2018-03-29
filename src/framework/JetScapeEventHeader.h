// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef JETSCAPEEVENTHEADER_H
#define JETSCAPEEVENTHEADER_H

/* #include "JetClass.h" */
/* #include "PartonShower.h" */

namespace Jetscape {

  /**
     Container for a multitude of event-related information
     such as xsec, centrality, ...
   */
  class JetScapeEventHeader {
    
  public:
    
  JetScapeEventHeader()
    :SigmaGen(0) {};
    // ~JetScapeEventHeader(){};
    // JetScapeEventHeader(const JetScapeEventHeader &c); //copy constructor
    
    /* const Parton& getParton(int idx) const; */
    /* const vector<Parton>& getPartonCollection() const; */
    /* void addParton(Parton &p); */
    /* void addPartonShower(shared_ptr<PartonShower> ps); */
    /* void deleteParton(int idx); */

    double GetSigmaGen(){return SigmaGen;};
    void SetSigmaGen( double d){ SigmaGen=d; };
    
    double GetSigmaErr(){return SigmaErr;};
    void SetSigmaErr( double d){ SigmaErr=d; };

    double GetEventWeight(){return EventWeight;};
    void SetEventWeight( double d){ EventWeight=d; };


  private:
    double SigmaGen;
    double SigmaErr;
    double EventWeight;
    
  };

} // end namespace Jetscape

#endif // JETSCAPEEVENTHEADER_H
