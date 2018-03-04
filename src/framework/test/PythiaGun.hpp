/** Create a pythia collision at a specified point and return the two inital hard partons
    
    Based on PGun and JSPythia

    @group JetScape (modular/task) based framework
    @author Kolja Kauder
    @version Revision 0.1
    @date Jun 29, 2017
*/

#ifndef PYTHIAGUN_HPP
#define PYTHIAGUN_HPP

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
};

#endif  // PYTHIAGUN_HPP
