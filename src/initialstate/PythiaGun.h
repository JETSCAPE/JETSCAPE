/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

// Create a pythia collision at a specified point and return the two inital hard partons

#ifndef PYTHIAGUN_H
#define PYTHIAGUN_H

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

class PythiaGun : public HardProcess, public Pythia8::Pythia {

private:
  double pTHatMin;
  double pTHatMax;
  double eCM;
  double vir_factor;
  bool initial_virtuality_pT;
  double softMomentumCutoff;
  bool FSR_on;
  bool softQCD;

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<PythiaGun> reg;

public:
  /** standard ctor
      @param xmlDir: Note that the environment variable PYTHIA8DATA takes precedence! So don't use it.
      @param printBanner: Suppress starting blurb. Should be set to true in production, credit where it's due
  */
  PythiaGun(string xmlDir = "DONTUSETHIS", bool printBanner = false)
      : Pythia8::Pythia(xmlDir, printBanner), HardProcess() {
    SetId("UninitializedPythiaGun");
  }

  ~PythiaGun();

  void InitTask();
  void Exec();
  
  double pthat_bins[66][2] = {
    {2400, 2510}, {2200, 2400}, {2000, 2200}, {1900, 2000}, {1800, 1900}, {1700, 1800}, {1600, 1700}, {1500, 1600},
    {1400, 1500}, {1300, 1400}, {1200, 1300}, {1100, 1200}, {1000, 1100}, {900, 1000}, {800, 900}, {700, 800},
    {600, 700}, {550, 600}, {500, 550}, {450, 500}, {400, 450}, {350, 400}, {300, 350}, {290, 300}, {280, 290},
    {270, 280}, {260, 270}, {250, 260}, {240, 250}, {230, 240}, {220, 230}, {210, 220}, {200, 210}, {190, 200},
    {180, 190}, {170, 180}, {160, 170}, {150, 160}, {140, 150}, {130, 140}, {120, 130}, {110, 120}, {100, 110},
    {90, 100}, {80, 90}, {70, 80}, {60, 70}, {55, 60}, {50, 55}, {45, 50}, {40, 45}, {35, 40}, {30, 35}, {25, 30},
    {20, 25}, {17, 20}, {15, 17}, {13, 15}, {11, 13}, {9, 11}, {7, 9}, {5, 7}, {4, 5}, {3, 4}, {2, 3}, {1, 2}
};


  // Getters
  double GetpTHatMin() const { return pTHatMin; }
  double GetpTHatMax() const { return pTHatMax; }

  // Cross-section information in mb and event weight.
  double GetSigmaGen() { return info.sigmaGen(); };
  double GetSigmaErr() { return info.sigmaErr(); };
  double GetPtHat() { return info.pTHat(); };
  double GetEventWeight() { return info.weight(); };
};

#endif // PYTHIAGUN_H
