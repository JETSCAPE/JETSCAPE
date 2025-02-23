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
// This is a general basic class for hadronic afterburner

#ifndef AFTERBURNER_H
#define AFTERBURNER_H

#include "JetScapeModuleBase.h"
#include "SoftParticlization.h"
#include "sigslot.h"

namespace Jetscape {

/**
 * @class Afterburner
 * @brief Interface to hadronic afterburner.
 *
 * The Afterburner class gathers hadrons from soft particlization and 
 * fragmentation, and performs additional processing before passing them to 
 * the final event output.
 */
class Afterburner : public JetScapeModuleBase {
 public:
  /**
   * @brief Default constructor.
   */
  Afterburner() {
    VERBOSE(8);
    SetId("Afterburner");
  }

  /**
   * @brief Destructor.
   */
  ~Afterburner() {
    VERBOSE(8);
    disconnect_all();
  }

  /**
   * @brief Initialize the Afterburner module.
   *
   * This function ensures that required configuration parameters are loaded
   * and initializes necessary random number distributions.
   */
  virtual void Init();

  /**
   * @brief Execute the afterburner process.
   *
   * This function runs the Afterburner module, handling the final state
   * processing of hadrons.
   */
  virtual void Exec();

 protected:
  /**
   * @brief Gather all hadrons from soft particlization and fragmentation.
   *
   * This function collects hadrons from different sources and prepares
   * them for the hadronic afterburner.
   *
   * @return A vector of vectors containing shared pointers to hadrons.
   */
  std::vector<std::vector<std::shared_ptr<Hadron>>> GatherAfterburnerHadrons();

  /**
   * @brief Retrieve soft particlization hadrons.
   *
   * This function obtains hadrons generated through the soft particlization
   * process.
   *
   * @return A vector of vectors containing shared pointers to soft particlization hadrons.
   */
  std::vector<std::vector<std::shared_ptr<Hadron>>>
  GetSoftParticlizationHadrons();

  /**
   * @brief Retrieve fragmentation hadrons.
   *
   * This function fetches hadrons generated from Hybrid Hadronization.
   *
   * @return A vector containing shared pointers to fragmentation hadrons.
   */
  std::vector<std::shared_ptr<Hadron>> GetFragmentationHadrons();

  std::vector<std::vector<std::shared_ptr<Hadron>>> dummy; ///< Placeholder storage for hadrons.
  std::uniform_real_distribution<double> ZeroOneDistribution; ///< Uniform random number distribution for position smearing.
  std::shared_ptr<std::uniform_int_distribution<int>> rand_int_ptr_; ///< RNG for Kaon-L / Kaon-S switching to K0 / Anti-K0.
};

}  // end namespace Jetscape

#endif  // AFTERBURNER_H
