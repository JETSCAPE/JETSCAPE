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

/**
 * @file JetScapeTaskSupport.h
 * @brief Declaration of the JetScapeTaskSupport singleton class.
 *
 * This class keeps track of every created task in a thread-safe manner and
 * provides shared or task-specific resources (like random number generators)
 * necessary for reproducible simulations. It acts as a support utility for
 * JetScapeTask modules.
 */

#ifndef JETSCAPETASKSUPPORT_H
#define JETSCAPETASKSUPPORT_H

#include <atomic>
#include <iostream>
#include <memory>
#include <random>
#include <thread>

#include "FluidDynamics.h"
#include "HardProcess.h"
#include "InitialState.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "JetScapeWriter.h"

using std::atomic_int;

namespace Jetscape {

/**
 * @class JetScapeTaskSupport
 * @brief Singleton class that provides support utilities for JetScape tasks.
 *
 * This class manages task registration, task-specific or shared random number
 * generators, and other resources. It ensures thread-safe access and consistent
 * initialization across tasks. 
 *
 * Initially developed to provide reproducible random seeds to each task.
 * It uses a Mersenne Twister engine (`std::mt19937`) as the RNG.
 *
 * - Uses magic statics (thread-safe except on MSVC 2013).
 * - `make_unique` was avoided due to platform-specific issues.
 */
class JetScapeTaskSupport {
 public:
  /**
   * @brief Get the singleton instance of JetScapeTaskSupport.
   * @return Pointer to the singleton instance.
   */
  static JetScapeTaskSupport *Instance();

  /**
   * @brief Register a task and assign it a unique ID.
   * 
   * This function should be called once per task at construction time.
   * @return An integer task ID unique to the current session.
   */
  int RegisterTask();

  /**
   * @brief Read the random seed from the XML configuration.
   * 
   * Initializes the random seed used for RNG generation.
   */
  static void ReadSeedFromXML();

  /**
   * @brief Get a random number generator for a task.
   *
   * Returns a `std::shared_ptr<std::mt19937>` to a Mersenne Twister RNG.
   * Depending on configuration, each task may receive a unique RNG or share a 
   * global one.
   *
   * @param TaskId The ID of the task requesting the RNG.
   * @return Shared pointer to a Mersenne Twister RNG.
   */
  shared_ptr<std::mt19937> GetMt19937Generator(int TaskId);

  /**
   * @brief Get the base random seed used for RNG creation.
   * @return The current random seed as an unsigned integer.
   */
  static unsigned int GetRandomSeed() { return random_seed_; };

 protected:
  /**
   * @brief Flag to determine if each task should have its own RNG.
   */
  static bool one_generator_per_task_;

 private:
  /**
   * @brief Private constructor (singleton pattern).
   */
  JetScapeTaskSupport() : CurrentTaskNumber(0){};

  /**
   * @brief Singleton instance.
   */
  static JetScapeTaskSupport *m_pInstance;

  /**
   * @brief Counter to assign unique task IDs.
   */
  atomic_int CurrentTaskNumber;

  /**
   * @brief The global random seed value.
   */
  static unsigned int random_seed_;

  /**
   * @brief Initialization flag for RNG.
   */
  static bool initialized_;

  /**
   * @brief Shared RNG instance for all tasks (if applicable).
   */
  static shared_ptr<std::mt19937> one_for_all_;
};

}  // end namespace Jetscape

#endif
