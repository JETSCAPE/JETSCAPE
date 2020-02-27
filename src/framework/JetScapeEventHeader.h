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

#ifndef JETSCAPEEVENTHEADER_H
#define JETSCAPEEVENTHEADER_H

namespace Jetscape {

/**
     Container for a multitude of event-related information
     such as xsec, centrality, ...
   */
class JetScapeEventHeader {

public:
  JetScapeEventHeader(){};
  // ~JetScapeEventHeader(){};
  // JetScapeEventHeader(const JetScapeEventHeader &c); //copy constructor

  /* const Parton& getParton(int idx) const; */
  /* const vector<Parton>& getPartonCollection() const; */
  /* void addParton(Parton &p); */
  /* void addPartonShower(shared_ptr<PartonShower> ps); */
  /* void deleteParton(int idx); */

  // ============================ Initial Hard Process =================================
  /// Initial Hard Process: Get cross section
  /// Note: In most cases, this value becomes more precise as more events are created.
  /// It is recommended to use the last event's value
  double GetSigmaGen() { return SigmaGen; };
  /// Initial Hard Process: Set cross section
  void SetSigmaGen(double d) { SigmaGen = d; };

  /// Initial Hard Process: Get uncertainty on the cross section
  /// Note: In most cases, this value becomes more smaller as more events are created.
  /// It is recommended to use the last event's value
  double GetSigmaErr() { return SigmaErr; };
  /// Initial Hard Process: Set uncertainty on the cross section
  void SetSigmaErr(double d) { SigmaErr = d; };

  /// Initial Hard Process: Get additionally created weight (e.g. pythia.event().weight())
  double GetEventWeight() { return EventWeight; };
  /// Initial Hard Process: Set additionally created weight (e.g. pythia.event().weight())
  void SetEventWeight(double d) { EventWeight = d; };

  // ============================ Initial State =================================
  /// Initial State: Get number of participants
  double GetNpart() { return Npart; };
  /// Initial State: Get number of participants
  void SetNpart(double d) { Npart = d; };

  /// Initial State: Get number of binary collisions
  double GetNcoll() { return Ncoll; };
  /// Initial State: Get number of binary collisions
  void SetNcoll(double d) { Ncoll = d; };

  /// Initial State: Get total entropy
  double GetTotalEntropy() { return TotalEntropy; };
  /// Initial State: Get total entropy
  void SetTotalEntropy(double d) { TotalEntropy = d; };

  // ============================ Hydro =================================
  /// Hydro: Get (2nd order) event plane angle
  double GetEventPlaneAngle() { return EventPlaneAngle; };
  /// Hydro: Set (2nd order) event plane angle
  void SetEventPlaneAngle(double d) { EventPlaneAngle = d; };

private:
  // ============================ Initial Hard Process =================================
  double SigmaGen = -1;
  double SigmaErr = -1;
  double EventWeight = 1;

  // ============================ Initial State =================================
  double Npart = -1; // could be int, but using double to allow averaged values
  double Ncoll = -1; // could be int, but using double to allow averaged values
  double TotalEntropy = -1;

  // ============================ Hydro =================================
  double EventPlaneAngle = -999;

  // ============================ Other possible options =================================
  // IS:
  // double Eccentricity;
  // double ImpactParameter;

  // Hydro:
  // angles, eccentricities

  // Eloss:
  // Switching criteria

  // Potential for consistency checks:
  // TotalEntropy from Free Streaming and from Hydro

  // Potential file-wide parameters (should be implemented in a different class)

  // IS:
  // string NuclearDensity; // (Woods-Saxon?)
  // string SaturationModel; // (IP-SAT, MCKLN, ...)

  // Inital Hard Process:
  // ptHat_min, ptHat_max
  // generator name, npdf,

  // Hydro:
  // EOS

  // Free Streaming:
  // Name, version

  // Hadronization:
  // type (colorless/colored, reco, ...)
};

} // end namespace Jetscape

#endif // JETSCAPEEVENTHEADER_H
