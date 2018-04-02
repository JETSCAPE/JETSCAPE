/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#ifndef JETSCAPEEVENTHEADER_H
#define JETSCAPEEVENTHEADER_H

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
