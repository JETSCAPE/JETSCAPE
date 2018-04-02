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

/** Create a pythia collision at a specified point and return the two inital hard partons
    @group JetScape (modular/task) based framework
    @version Revision 0.1
    @date Jun 29, 2017
*/

#ifndef PYTHIAGUN_H
#define PYTHIAGUN_H

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

class PythiaGun: public HardProcess, public Pythia8::Pythia {
   
private:
  double pTHatMin;
  double pTHatMax;
  double eCM;

public:
  /** standard ctor
      @param xmlDir: Note that the environment variable PYTHIA8DATA takes precedence! So don't use it.
      @param printBanner: Suppress starting blurb. Should be set to true in production, credit where it's due
  */
  PythiaGun(string xmlDir = "DONTUSETHIS", bool printBanner = false)
    : Pythia8::Pythia(xmlDir,printBanner), HardProcess()
  {
    SetId("UninitializedPythiaGun");
  }
  
  ~PythiaGun();
 
  void InitTask();
  void Exec();

  // Getters
  double GetpTHatMin() const {return pTHatMin;}
  double GetpTHatMax() const {return pTHatMax;}

  // Cross-section information in mb and event weight.
  double GetSigmaGen(){ return info.sigmaGen();};
  double GetSigmaErr(){ return info.sigmaErr();};  
  double GetEventWeight(){ return info.weight(); };
  
};

#endif  // PYTHIAGUN_H
