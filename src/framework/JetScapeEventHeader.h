/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
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
 * @class JetScapeEventHeader
 * @brief Container for event-related information such as cross-section, centrality, etc.
 */
 class JetScapeEventHeader {
 public:

  /**
   * @brief Default constructor
   */
  JetScapeEventHeader(){};

  // ~JetScapeEventHeader(){};
  // JetScapeEventHeader(const JetScapeEventHeader &c); //copy constructor

  /* const Parton& getParton(int idx) const; */
  /* const vector<Parton>& getPartonCollection() const; */
  /* void addParton(Parton &p); */
  /* void addPartonShower(shared_ptr<PartonShower> ps); */
  /* void deleteParton(int idx); */

  // ============================ Initial Hard Process
  // =================================

  /**
   * @brief Initial Hard Process: Get cross-section.
   * @note In most cases, this value becomes more precise as more events are
   * created. It is recommended to use the last event's value.
   * @return Cross-section value
   */
  double GetSigmaGen() { return SigmaGen; };

   /**
   * @brief Initial Hard Process: Set the generated cross-section
   * @param d Cross-section value
   */
  void SetSigmaGen(double d) { SigmaGen = d; };

  /**
   * @brief Initial Hard Process: Get uncertainty on the cross-section
   * @note Note: In most cases, this value becomes more smaller as more
   * events are created. It is recommended to use the last event's value.
   * @return Cross-section uncertainty value
   */
  double GetSigmaErr() { return SigmaErr; };

  /**
   * @brief Initial Hard Process: Set uncertainty on the cross-section
   * @param d Cross-section uncertainty value
   */
  void SetSigmaErr(double d) { SigmaErr = d; };

  /**
   * @brief Initial Hard Process: Get the pt-hat value
   * @return pt-hat value
   */
  double GetPtHat() { return PtHat; };

  /**
   * @brief Set the pt-hat value
   * @param d pt-hat value
   */
  void SetPtHat(double d) { PtHat = d; };

  /**
   * @brief Initial Hard Process: Get additionally created weight (e.g. pythia.event().weight())
   * @return Event weight
   */
  double GetEventWeight() { return EventWeight; };

   /**
   * @brief Initial Hard Process: Set additionally created weight (e.g. pythia.event().weight())
   * @param d Event weight
   */
  void SetEventWeight(double d) { EventWeight = d; };

  // ============================ Initial State
  // =================================

  /**
   * @brief Initial State: Get number of participants
   * @return Number of participants
   */
  double GetNpart() { return Npart; };

  /**
   * @brief Initial State: Set number of participants
   * @param d Number of participants
   */
  void SetNpart(double d) { Npart = d; };

  /**
   * @brief Initial State: Get the number of binary collisions
   * @return Number of binary collisions
   */
  double GetNcoll() { return Ncoll; };

  /**
   * @brief Initial State: Set the number of binary collisions
   * @param d Number of binary collisions
   */
  void SetNcoll(double d) { Ncoll = d; };

  /**
   * @brief Initial State: Get centrality of the event
   * @return Event centrality
   */
  double GetEventCentrality() { return EventCentrality; }

  /**
   * @brief Initial State: Set the event centrality
   * @param d Event centrality
   */
  void SetEventCentrality(double d) { EventCentrality = d; }

  /**
   * @brief  Initial State: Get the total entropy
   * @return Total entropy
   */
  double GetTotalEntropy() { return TotalEntropy; };

   /**
   * @brief Initial State: Set the total entropy
   * @param d Total entropy
   */
  void SetTotalEntropy(double d) { TotalEntropy = d; };

  // ============================ Hydro =================================

  /**
   * @brief Hydro: Get the event plane angle (2nd order)
   * @return Event plane angle
   */
  double GetEventPlaneAngle() { return EventPlaneAngle; };

  /**
   * @brief Hydro: Set the event plane angle (2nd order)
   * @param d Event plane angle
   */
  void SetEventPlaneAngle(double d) { EventPlaneAngle = d; };

 private:

  // ============================ Initial Hard Process
  // =================================

  double SigmaGen = -1;          ///< Cross-section value
  double SigmaErr = -1;          ///< Cross-section uncertainty
  double PtHat = -1;             ///< pt-hat value
  double EventWeight = 1;        ///< Event weight

  // ============================ Initial State
  // =================================

  double Npart = -1;             ///< Number of participants (could be int, but, using double to allow averaged values)
  double Ncoll = -1;             ///< Number of binary collisions (could be int, but using double to allow averaged values)
  double TotalEntropy = -1;      ///< Total entropy
  double EventCentrality = -1;   ///< Event centrality

  // ============================ Hydro =================================

  double EventPlaneAngle = -999; ///< Event plane angle (2nd order)

  // ============================ Other possible options
  // ================================= IS: double Eccentricity; double
  // ImpactParameter;

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

}  // end namespace Jetscape

#endif  // JETSCAPEEVENTHEADER_H
