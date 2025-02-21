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

#ifndef HADRONIZATIONMANAGER_H
#define HADRONIZATIONMANAGER_H

#include <vector>

#include "Hadronization.h"
#include "JetClass.h"
#include "JetScapeTask.h"
#include "sigslot.h"

namespace Jetscape {

/**
 * @class HadronizationManager
 * @brief Manages hadronization tasks within the JETSCAPE framework.
 *
 * This class is responsible for coordinating hadronization submodules,
 * handling signals, and managing the transformation of partons into hadrons.
 */
class HadronizationManager
    : public JetScapeModuleBase,
      public std::enable_shared_from_this<HadronizationManager> {
 public:
  /**
   * @brief Default constructor.
   */
  HadronizationManager();
  
  /**
   * @brief Destructor.
   */
  virtual ~HadronizationManager();

  /**
   * @brief Initializes the hadronization manager.
   */
  virtual void Init();

  /**
   * @brief Executes the hadronization process.
   */
  virtual void Exec();

  /**
   * @brief Clears internal states and resources.
   */
  virtual void Clear();

  /**
   * @brief Writes task-specific data to an output writer.
   * @param w Weak pointer to a JetScapeWriter.
   */
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);

  /**
   * @brief Retrieves the number of registered signals.
   * @return Number of signals.
   */
  int GetNumSignals();

  /**
   * @brief Creates the required signal slots for hadronization modules.
   */
  void CreateSignalSlots();

  /**
   * @brief Retrieves hadrons from hadronization submodules.
   * @param signal Reference to a vector where hadrons will be stored.
   */
  void GetHadrons(vector<shared_ptr<Hadron>> &signal);

  /**
   * @brief Deletes hadrons from hadronization submodules.
   * 
   * This function is used when hadrons are passed to the afterburner.
   * Otherwise, hadrons are printed to file and the same hadrons are modified
   * in the transport and printed again.
   */
  void DeleteHadrons();

  /**
   * @brief Deletes hadrons with a positive status flag.
   * 
   * This function is used when hadrons are passed to the afterburner. The 
   * negative status flag hadrons are not deleted (can not be propagated in the
   * afterburner).
   */
  void DeleteRealHadrons();

  /**
   * Signal to retrieve hadrons from the hard process, not from hadronization submodules.
   */
  sigslot::signal1<vector<shared_ptr<Hadron>> &>
      GetHadronList;
  
  /**
   * Signal to retrieve the final list of partons before hadronization.
   */
  sigslot::signal1<vector<vector<shared_ptr<Parton>>> &> GetFinalPartonList;

  /**
   * @brief Sets the connection status of the final parton list signal.
   * @param m_GetFinalPartonListConnected Connection status (true/false).
   */
  void SetGetFinalPartonListConnected(bool m_GetFinalPartonListConnected) {
    GetFinalPartonListConnected = m_GetFinalPartonListConnected;
  }

  /**
   * @brief Gets the connection status of the final parton list signal.
   * @return True if connected, false otherwise.
   */
  const bool GetGetFinalPartonListConnected() {
    return GetFinalPartonListConnected;
  }

  /**
   * @brief Sets the connection status of the hadron list signal.
   * @param m_GetHadronListConnected Connection status (true/false).
   */
  void SetGetHadronListConnected(bool m_GetHadronListConnected) {
    GetHadronListConnected = m_GetHadronListConnected;
  }

  /**
   * @brief Gets the connection status of the hadron list signal.
   * @return True if connected, false otherwise.
   */
  const bool GetGetHadronListConnected() { return GetHadronListConnected; }

 private:
  /// Connection status of the final parton list signal.
  bool GetFinalPartonListConnected;

  /// Connection status of the hadron list signal.
  bool GetHadronListConnected;

  /// Final parton list before hadronization.
  vector<vector<shared_ptr<Parton>>> hd;

  /// Hadron list from hadronization submodules.
  vector<shared_ptr<Hadron>> hadrons;
};

}  // end namespace Jetscape

#endif
