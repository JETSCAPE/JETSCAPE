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

#ifndef EPEMGUN_H
#define EPEMGUN_H

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

using namespace Jetscape;

class epemGun : public HardProcess, public Pythia8::Pythia {

private:
  //double pTHatMin;
  //double pTHatMax;
  double eCM;
  //bool FSR_on;

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<epemGun> reg;

  double sud_val_QG_w_M(double M, double h0, double h1, double h2, double E1);
  double sud_z_QG_w_M(double M, double cg, double cg1, double E2);
  double alpha_s(double q2);

protected:
  std::uniform_real_distribution<double> ZeroOneDistribution;

public:
  /** standard ctor
      @param xmlDir: Note that the environment variable PYTHIA8DATA takes precedence! So don't use it.
      @param printBanner: Suppress starting blurb. Should be set to true in production, credit where it's due
  */
  epemGun(string xmlDir = "DONTUSETHIS", bool printBanner = false)
      : Pythia8::Pythia(xmlDir, printBanner), HardProcess() {
    SetId("UninitializedepemGun");
  }

  ~epemGun();

  void InitTask();
  void Exec();

  // Getters
  //double GetpTHatMin() const { return pTHatMin; }
  //double GetpTHatMax() const { return pTHatMax; }

  // Cross-section information in mb and event weight.
  double GetSigmaGen() { return info.sigmaGen(); };
  double GetSigmaErr() { return info.sigmaErr(); };
  double GetEventWeight() { return info.weight(); };
};

#endif // EPEMGUN_H
